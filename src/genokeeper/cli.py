"""The `genokeeper` command line.

Each sub-command lives in its own module and joins the parser through a
`register(subparsers)` function listed in `REGISTRARS`. A registrar adds its sub-parser and
sets a `handler` default: a function from the parsed arguments to an exit code.
"""

import argparse
import sys
from collections.abc import Callable

from genokeeper import (
    __version__,
    carrier,
    changes,
    check,
    dashboard,
    explain,
    pipeline,
    query,
    refresh,
    schedule,
)
from genokeeper.calibrate import build as calibrate
from genokeeper.evidence import frequency

type Subparsers = argparse._SubParsersAction[argparse.ArgumentParser]
type Registrar = Callable[[Subparsers], None]
type Handler = Callable[[argparse.Namespace], int]

REGISTRARS: tuple[Registrar, ...] = (
    query.register,
    carrier.register,
    changes.register,
    check.register,
    dashboard.register,
    explain.register,
    pipeline.register,
    refresh.register,
    schedule.register,
    calibrate.register,
    frequency.register,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genokeeper",
        description="Findings from your own whole-genome data, kept on this machine.",
    )
    parser.add_argument(
        "--version", action="store_true", help="print the version and exit"
    )
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    for register in REGISTRARS:
        register(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.version:
        print(f"genokeeper {__version__}")
        return 0
    handler: Handler | None = getattr(args, "handler", None)
    if handler is None:
        parser.print_usage(sys.stderr)
        return 2
    return handler(args)
