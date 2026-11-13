"""`genokeeper demo [--dir PATH]`: download the packaged HG002 store from the demo's release
(`DEMO_BASE`) and write a folder Claude Code can work in.

`MANIFEST` comes first and is checked before any store byte: its schema version must be this
genokeeper's `SCHEMA_VERSION`. Each file is streamed into `<name>.part` with its sha256
computed on the way and renamed to its final name only when that matches the manifest's. The
sha256 catches a cut or corrupted download; the anchor for what is downloaded is HTTPS to
GitHub's release, which a manifest from the same origin cannot add to.

The folder (`install`'s docstring lists it) is marked by `DEMO_MARKER`, so `load_paths` reads
its own `lab/` whatever `$GENOKEEPER_HOME` says. An existing folder whose marker is not a demo
manifest is refused; one with it is updated in place, and files of anyone else's in it stay.
"""

import argparse
import hashlib
import http.client
import json
import os
import shlex
import shutil
import sys
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from genokeeper import __version__, refresh_log
from genokeeper.config import DEMO_MARKER, load_paths, skill_resource
from genokeeper.dashboard import write_dashboard
from genokeeper.demo import MANIFEST, REFRESH_LOG_FILE, STORE_FILE
from genokeeper.demo.frame import CALIBRATION_NOTE, FRAME, FRAME_SUMMARY
from genokeeper.query import RunNotFound, open_store
from genokeeper.store.db import SCHEMA_VERSION, connect, export_snapshots

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

DEMO_BASE = "https://github.com/genokeeper/hg002-demo/releases/download/hg002/"
TIMEOUT = 60  # seconds for any one network read
_CHUNK = 1 << 20
_MANIFEST_MAX = 1 << 20
_USER_AGENT = f"genokeeper/{__version__} (demo download)"
_SAY = "genokeeper demo"
FOLDER = "~/genokeeper-demo"  # the folder `genokeeper demo` writes by default

QUESTIONS = (
    "Is HG002 a carrier of anything their children could inherit?",
    "How would HG002 respond to clopidogrel, which CYP2C19 activates?",
    "List every change the last ClinVar update brought, below the default thresholds too.",
)
# What Finder writes into a folder it opens: such a folder is not "something else".
_FINDER = (".DS_Store", ".localized")
# What a question cannot get here, each with its reason; CLAUDE.md lists them.
CANNOT = (
    (
        "Plain-language explanations are not generated here: `genokeeper explain` writes"
        " them with Claude from MedlinePlus files this folder does not have, so the"
        " dashboard shows none."
    ),
)

_PYPROJECT = """\
[project]
name = "genokeeper-demo"
version = "0"
requires-python = ">=3.14"
dependencies = ["genokeeper=={version}"]
"""


class DemoRefused(Exception):
    """The demo was not written: the folder is someone else's, the download failed or does
    not match its manifest, or the store's schema is not this genokeeper's."""


def _open(url: str) -> http.client.HTTPResponse:
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    return urllib.request.urlopen(request, timeout=TIMEOUT)


def _manifest(base: str) -> tuple[bytes, dict[str, str]]:
    """The manifest's bytes and the sha256 it names for the store and the refresh log;
    refuse one whose schema version is not this genokeeper's."""
    url = base + MANIFEST
    try:
        with _open(url) as response:
            body = response.read(_MANIFEST_MAX)
    except (OSError, http.client.HTTPException) as e:
        raise DemoRefused(f"could not download {url}: {e}") from e
    try:
        manifest: dict[str, Any] = json.loads(body)
        schema = manifest["schema_version"]
        sums = {
            name: str(manifest["files"][name])
            for name in (STORE_FILE, REFRESH_LOG_FILE)
        }
    except (ValueError, KeyError, TypeError) as e:
        raise DemoRefused(f"{url} is not a demo manifest: {e!r}") from e
    if schema != SCHEMA_VERSION:
        made = manifest.get("version", "another version")
        how = (
            "upgrade genokeeper (uvx genokeeper@latest demo)"
            if isinstance(schema, int) and schema > SCHEMA_VERSION
            else f"run the genokeeper that packaged it (uvx genokeeper=={made} demo)"
        )
        raise DemoRefused(
            f"the demo's store has schema version {schema}, and this genokeeper"
            f" ({__version__}) reads {SCHEMA_VERSION}: {how}"
        )
    return body, sums


_VIEWS = "SELECT view_name, sql FROM duckdb_views() WHERE NOT internal ORDER BY 1"
_MACROS = "SELECT function_name FROM duckdb_functions() WHERE NOT internal"
_INERT: dict[str, str | bool | int | float | list[str]] = {
    "enable_external_access": False,
    "autoload_known_extensions": False,
    "autoinstall_known_extensions": False,
}


def _inert(store: Path) -> None:
    """Refuse (`DemoRefused`) a store whose views are not exactly the ones this genokeeper's
    schema creates, or that holds a macro: a view runs SQL on every later read, and could read
    the user's files. The adversary is whoever can replace the release's assets (the manifest
    comes from the same release, so the sha256 cannot stop it); the anchor is the expected
    views, from the installed wheel's schema, which the asset cannot move, and the catalog read
    with external access off, which runs no view."""
    try:
        with duckdb.connect(str(store), read_only=True, config=_INERT) as con:
            views, macros = (
                con.execute(_VIEWS).fetchall(),
                con.execute(_MACROS).fetchall(),
            )
    except duckdb.Error as e:
        raise DemoRefused(f"the demo's store cannot be read: {e}") from e
    with connect(Path(":memory:")) as fresh:
        expected = fresh.execute(_VIEWS).fetchall()
    if views != expected or macros:
        raise DemoRefused(
            "the demo's store holds views or macros this genokeeper's schema does not"
            " create; nothing was installed"
        )


def _download(
    url: str,
    sha256: str,
    target: Path,
    *,
    size_line: bool = False,
    check: Callable[[Path], None] | None = None,
) -> None:
    """`url` into `target`, through `<target>.part`, renamed only when its sha256 is
    `sha256` and `check`, given the part, raises nothing. With `size_line`, say its size
    first, and show progress on a terminal."""
    part = target.with_name(target.name + ".part")
    target.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    progress = size_line and sys.stderr.isatty()
    try:
        try:
            with _open(url) as response, part.open("wb") as out:
                length = response.headers.get("Content-Length")
                total = int(length) if length and length.isdigit() else None
                if size_line:
                    size = "size unknown" if total is None else f"{total / 1e6:.1f} MB"
                    print(f"{_SAY}: downloading HG002's store, {size}", file=sys.stderr)
                done = 0
                while chunk := response.read(_CHUNK):
                    out.write(chunk)
                    digest.update(chunk)
                    done += len(chunk)
                    if progress:
                        of = "" if total is None else f" of {total / 1e6:.1f}"
                        print(f"\r  {done / 1e6:.1f}{of} MB", end="", file=sys.stderr)
            if progress:
                print(file=sys.stderr)
        except (OSError, http.client.HTTPException) as e:
            raise DemoRefused(f"could not download {url}: {e}") from e
        if digest.hexdigest() != sha256:
            raise DemoRefused(
                f"{url}: its sha256 does not match the manifest's (a cut or corrupted"
                " download); run genokeeper demo again"
            )
        if check is not None:
            check(part)
        os.replace(part, target)
    finally:
        part.unlink(missing_ok=True)


def _claude_md() -> str:
    cannot = "\n".join(f"- {line}" for line in CANNOT)
    questions = "\n".join(f"{n}. {q}" for n, q in enumerate(QUESTIONS, 1))
    return f"""\
# HG002 demo

`genokeeper demo` rewrites this file; keep your own notes in `CLAUDE.local.md`.

{FRAME_SUMMARY}

{FRAME}

{CALIBRATION_NOTE}

In your first answer here, say once: "{FRAME_SUMMARY}"

The genome in this folder is HG002's. Answer questions about it ("HG002", "this genome",
"my genome") with the genome-ask skill (`.claude/skills/genome-ask/SKILL.md`): its
`uv run genokeeper …` commands read this folder's store (`reports/genomelab.duckdb`) and its
own lab (`lab/`), with the genokeeper version `pyproject.toml` pins.

## What cannot be answered here

{cannot}

## Questions to try

{questions}
"""


def _ours(folder: Path) -> bool:
    """`folder`'s `DEMO_MARKER` is a demo manifest (it has `schema_version` and `runs`)."""
    try:
        marker = json.loads((folder / DEMO_MARKER).read_bytes())
    except OSError, ValueError:
        return False
    return isinstance(marker, dict) and {"schema_version", "runs"} <= marker.keys()


def _foreign(folder: Path) -> bool:
    """`folder` holds anything but Finder's files and empty folders (a refused first
    download leaves `reports/`)."""
    return any(
        p.name not in _FINDER
        and not p.name.startswith("._")
        and not (p.is_dir() and not any(p.iterdir()))
        for p in folder.iterdir()
    )


def install(folder: Path, *, base: str = DEMO_BASE) -> Path:
    """Write the demo folder at `folder` from the release at `base`; return `folder`:
    `DEMO_MARKER` (the manifest), `reports/` (the store and its snapshots), `lab/` (the
    marker and the refresh-log record), `dashboard.html` with the frame, `pyproject.toml`
    (this genokeeper, for the skill's `uv run`), `.claude/skills/genome-ask/SKILL.md` and
    `CLAUDE.md`. Raises `DemoRefused`; nothing is written before the manifest's check, and
    no file is under its final name unless its sha256 matched."""
    holds = folder.is_dir() and _foreign(folder)
    if (folder.exists() and not folder.is_dir()) or (holds and not _ours(folder)):
        raise DemoRefused(
            f"{folder} holds something else and is not a demo folder (no {DEMO_MARKER});"
            " choose another with --dir"
        )
    manifest, sums = _manifest(base)
    folder.mkdir(parents=True, exist_ok=True)
    store = folder / "reports" / STORE_FILE
    _download(base + STORE_FILE, sums[STORE_FILE], store, size_line=True, check=_inert)
    (folder / DEMO_MARKER).write_bytes(manifest)
    paths = load_paths(folder)  # the marker makes the folder's own `lab/` the lab
    log = paths.lab / "refresh" / refresh_log.LOG
    _download(base + REFRESH_LOG_FILE, sums[REFRESH_LOG_FILE], log)
    (paths.lab / DEMO_MARKER).write_text("{}\n", encoding="utf-8")
    snapshots = paths.reports / "snapshots"
    shutil.rmtree(snapshots, ignore_errors=True)
    print(f"{_SAY}: writing the store's snapshots", file=sys.stderr)
    with open_store(store) as con:
        export_snapshots(con, snapshots)
    print(f"{_SAY}: writing the dashboard, a few minutes", file=sys.stderr)
    write_dashboard(
        paths,
        out=folder / "dashboard.html",
        lang="en",
        banner=(FRAME_SUMMARY, FRAME, CALIBRATION_NOTE),
    )
    (folder / "pyproject.toml").write_text(
        _PYPROJECT.format(version=__version__), encoding="utf-8"
    )
    skill = folder / ".claude" / "skills" / "genome-ask" / "SKILL.md"
    skill.parent.mkdir(parents=True, exist_ok=True)
    skill.write_bytes(skill_resource().read_bytes())
    (folder / "CLAUDE.md").write_text(_claude_md(), encoding="utf-8")
    return folder


def _command(args: argparse.Namespace) -> int:
    folder = Path(args.dir).expanduser().absolute()
    try:
        install(folder, base=DEMO_BASE)
    except (DemoRefused, RunNotFound, OSError) as e:
        print(f"{_SAY}: {e}", file=sys.stderr)
        return 2
    print(f"{_SAY}: wrote {folder}")
    print(f"The dashboard, to open in a browser: {folder / 'dashboard.html'}")
    print("To ask Claude Code about HG002:")
    print(f"  cd {shlex.quote(str(folder))} && claude")
    print("Questions to try:")
    for n, question in enumerate(QUESTIONS, 1):
        print(f"  {n}. {question}")
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `demo` sub-command."""
    d = subparsers.add_parser(
        "demo", help="download HG002's findings into a folder to ask Claude Code about"
    )
    d.add_argument(
        "--dir",
        type=Path,
        default=Path(FOLDER),
        help=f"the folder to write (default: {FOLDER})",
    )
    d.set_defaults(handler=_command)
