"""The analyser's identity: what a run's findings depend on besides the user's data and
the knowledge-base releases.

`analyser_identity` is a sha256 over the package's code and resources (every file under
`PACKAGE_ROOT`, bytecode apart) and config/grading.toml. It is recorded in every run as
`run.tools["analyser"]`; a refresh compares it with the latest snapshot's to tell a change
of the analyser from a change of the sources. It is in no step key: the cached slow steps
are keyed by their own inputs, so a new analyser re-grades without re-calling.
"""

import hashlib
from importlib import resources
from pathlib import Path

from genokeeper.config import GRADING_FILE

PACKAGE_ROOT = Path(str(resources.files("genokeeper")))


def _files() -> list[tuple[str, Path]]:
    """Every file to hash as (name, path), in sorted name order: the package's by path
    relative to its root, then the grading config's."""
    package = sorted(
        (path.relative_to(PACKAGE_ROOT).as_posix(), path)
        for path in PACKAGE_ROOT.rglob("*")
        if path.is_file()
        and path.suffix != ".pyc"
        and "__pycache__" not in path.relative_to(PACKAGE_ROOT).parts
        and path.name != ".DS_Store"
        # The skill is documentation for the reader; editing it is not a new analyser.
        and path.relative_to(PACKAGE_ROOT).parts[:2] != ("resources", "skills")
    )
    return [*package, (f"config/{GRADING_FILE.name}", GRADING_FILE)]


def analyser_identity() -> str:
    """The hex sha256 of each file as `name\\0length\\0bytes`."""
    digest = hashlib.sha256()
    for name, path in _files():
        data = path.read_bytes()
        digest.update(f"{name}\0{len(data)}\0".encode() + data)
    return digest.hexdigest()
