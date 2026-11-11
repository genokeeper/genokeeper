"""`python -m genokeeper.demo package OUT`: the maintainer's demo commands, kept out of the
users' CLI (see `genokeeper.demo`)."""

import argparse
import sys
from pathlib import Path

from genokeeper.config import load_paths
from genokeeper.demo.package import Leak, package, private_id_file

# Where the files go once written: the demo repository's moving release, then the Pages
# workflow that checks them against the manifest and deploys them.
_PUBLISH = (
    "gh release upload hg002 {out}/* --clobber --repo genokeeper/hg002-demo",
    "gh workflow run pages.yml --repo genokeeper/genokeeper",
)


def _package(args: argparse.Namespace) -> int:
    try:
        out = package(
            load_paths().lab,
            args.out,
            home=Path.home(),
            private_ids=private_id_file(),
        )
    except (Leak, OSError) as e:
        print(f"genokeeper demo package: {e}", file=sys.stderr)
        return 2
    print(f"genokeeper demo package: wrote {out}; to publish:")
    for command in _PUBLISH:
        print(f"  {command.format(out=out)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m genokeeper.demo")
    commands = parser.add_subparsers(required=True)
    p = commands.add_parser(
        "package", help="write the demo's store, pages and manifest"
    )
    p.add_argument("out", type=Path, help="the folder to write into")
    p.set_defaults(handler=_package)
    args = parser.parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    sys.exit(main())
