"""`python -m genokeeper.demo build` and `package OUT`: the maintainer's demo commands, kept
out of the users' CLI (see `genokeeper.demo`)."""

import argparse
import shlex
import sys
from pathlib import Path

from genokeeper import __version__
from genokeeper.checksum import sha256_file
from genokeeper.config import load_paths
from genokeeper.demo import MANIFEST, SHIPPED
from genokeeper.demo.build import BuildRefused, build
from genokeeper.demo.package import Leak, package, private_id_file
from genokeeper.pipeline import RunRefused

# Where the files go once written, each named (never OUT/*): a release of their own, never
# replaced (an installed genokeeper pins it), its pin, then the Pages workflow that checks them
# against it and deploys them.
_PUBLISH = (
    (
        "publish",
        "gh release create {tag} --repo genokeeper/hg002-demo --title {tag} {files}",
    ),
    ("pin, committed before the wheel is built", "echo {pin} > {pin_file}"),
    ("then, from main", "gh workflow run pages.yml --repo genokeeper/genokeeper"),
)
_RELEASES = "https://github.com/genokeeper/hg002-demo/releases/download/"


def _package(args: argparse.Namespace) -> int:
    try:
        out = package(
            load_paths().lab,
            args.out,
            home=Path.home(),
            private_ids=private_id_file(),
            no_private_ids=args.no_private_ids,
        )
    except (Leak, OSError) as e:
        print(f"genokeeper demo package: {e}", file=sys.stderr)
        return 2
    print(f"genokeeper demo package: wrote {out}; to publish:")
    files = " ".join(shlex.quote(str(out / name)) for name in (MANIFEST, *SHIPPED))
    tag = f"hg002-{__version__}"
    pin = shlex.quote(f"{_RELEASES}{tag}/ {sha256_file(out / MANIFEST)}")
    pin_file = "src/genokeeper/resources/demo/release.pin"
    for when, command in _PUBLISH:
        said = command.format(files=files, tag=tag, pin=pin, pin_file=pin_file)
        print(f"  {when}: {said}")
    return 0


def _build(args: argparse.Namespace) -> int:
    try:
        first, second = build(load_paths().lab, threads=args.threads)
    except (BuildRefused, RunRefused, OSError) as e:
        print(f"genokeeper demo build: {e}", file=sys.stderr)
        return 2
    print(f"genokeeper demo build: runs {first} and {second}; next, `package OUT`")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m genokeeper.demo")
    commands = parser.add_subparsers(required=True)
    p = commands.add_parser(
        "package", help="write the demo's store, pages and manifest"
    )
    p.add_argument("out", type=Path, help="the folder to write into")
    p.add_argument(
        "--no-private-ids",
        action="store_true",
        help="package without the private-ID file's check (refused without it)",
    )
    p.set_defaults(handler=_package)
    b = commands.add_parser(
        "build", help="run HG002 on two ClinVar releases into the demo root"
    )
    b.add_argument("--threads", type=int, default=None, help="CPU threads")
    b.set_defaults(handler=_build)
    args = parser.parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    sys.exit(main())
