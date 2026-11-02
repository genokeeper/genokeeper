"""The daily update check as a user LaunchAgent, and the notification it sends (spec, "The
daily check").

`genokeeper schedule install` writes `~/Library/LaunchAgents/local.genokeeper.check-updates.plist`,
loads it (`launchctl bootout`, `enable`, then `bootstrap gui/<uid>`: bootstrapping a loaded
or a disabled label fails) and sends one test notification, whose exit status decides whether the install
worked. A user agent runs with `PATH=/usr/bin:/bin:/usr/sbin:/sbin`, so the plist names
absolute paths (`uv`, the notifier) and the main checkout as the project, never a worktree
that may be removed, and the install's language (`genokeeper.lang`), which the notifications
are written in. Background Task Management can still block an agent that launchd calls
enabled, so nothing trusts launchd's status: `check.json`'s timestamp says whether the check
runs (`STALE_DAYS`).

`notify` hands the text to terminal-notifier as an argument (`-message <text>`); upstream
strings are untrusted, and the text is never part of a script. terminal-notifier is required:
it reports a failure by its exit code, while `osascript display notification` can exit 0 and
show nothing, so without it the install could prove nothing and is refused.
"""

import argparse
import os
import plistlib
import shutil
import subprocess
import sys
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genokeeper.config import load_paths
from genokeeper.env import getenv
from genokeeper.lang import ENV as LANG_ENV
from genokeeper.lang import Language, language

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

LABEL = "local.genokeeper.check-updates"
# The agent's label before the rename from genomelab: its plist runs `genomelab
# check-updates`, which no longer exists, so `install` and `uninstall` remove it too.
OLD_LABEL = "local.genomelab.check-updates"
STALE_DAYS = 35  # a check that has not run for longer than this is silenced
# The absolute path of terminal-notifier, from install (`$GENOMELAB_NOTIFIER` before the
# rename, which an agent installed before it still sets).
NOTIFIER_ENV = "GENOKEEPER_NOTIFIER"
TITLE = "Genokeeper"
TEST_TEXT: dict[Language, str] = {
    "en": "Update checks are on: this is a test notification.",
    "ru": "Проверка обновлений включена: это тестовое уведомление.",
}
NO_NOTIFIER = 127  # `notify`'s status without terminal-notifier: "command not found"

type Launchctl = Callable[[list[str]], int]
type Notifier = Callable[[str], int]


class ScheduleError(RuntimeError):
    """The agent could not be installed."""


def main_checkout(repo: Path) -> Path:
    """The main checkout of the repository `repo` belongs to (`repo` may be a worktree,
    which can be removed): the parent of git's common directory."""
    done = subprocess.run(
        [
            "git",
            "-C",
            str(repo),
            "rev-parse",
            "--path-format=absolute",
            "--git-common-dir",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    common = Path(done.stdout.strip())
    if done.returncode != 0 or common.name != ".git":
        raise ScheduleError(f"{repo} is not in a git checkout with a .git directory")
    return common.parent


def agent_path(home: Path, label: str = LABEL) -> Path:
    return home / "Library" / "LaunchAgents" / f"{label}.plist"


def plist(
    repo: Path,
    lab: Path,
    uv: Path,
    *,
    terminal_notifier: Path,
    lang: Language = "en",
) -> dict[str, Any]:
    """The agent: `uv run --project <main checkout> --frozen genokeeper check-updates
    --notify` every day at 10:00 and at load, absolute paths throughout, its notifications
    in `lang` and its output in `<lab>/refresh/check.log`."""
    log = str(lab / "refresh" / "check.log")
    environment = {
        "GENOKEEPER_HOME": str(lab),
        NOTIFIER_ENV: str(terminal_notifier),
        LANG_ENV: lang,
    }
    return {
        "Label": LABEL,
        "ProgramArguments": [
            str(uv),
            "run",
            "--project",
            str(main_checkout(repo)),
            "--frozen",
            "genokeeper",
            "check-updates",
            "--notify",
        ],
        "StartCalendarInterval": {"Hour": 10, "Minute": 0},
        "RunAtLoad": True,
        "EnvironmentVariables": environment,
        "StandardOutPath": log,
        "StandardErrorPath": log,
    }


def notify(
    text: str,
    *,
    terminal_notifier: Path | None = None,
    environ: Mapping[str, str] = os.environ,
    run: Callable[..., subprocess.CompletedProcess[bytes]] = subprocess.run,
) -> int:
    """Show `text` as a notification and return the notifier's exit status (`NO_NOTIFIER`
    without one). The text is one argument of `terminal-notifier` (the path `install` wrote
    into the plist's environment), never part of a script."""
    found = terminal_notifier or (
        Path(path) if (path := getenv("NOTIFIER", environ)) else None
    )
    if found is None:
        return NO_NOTIFIER
    argv = [str(found), "-title", TITLE, "-message", text]
    done = run(argv, check=False, capture_output=True, stdin=subprocess.DEVNULL)
    if done.returncode != 0 and done.stderr:
        # The notifier's own reason, e.g. "Notifications are not allowed for this application".
        sys.stderr.write(done.stderr.decode(errors="replace"))
    return done.returncode


def launchctl(args: list[str]) -> int:
    """`launchctl <args>`, its output discarded: only the exit status is read."""
    done = subprocess.run(
        ["/bin/launchctl", *args],
        check=False,
        capture_output=True,
        stdin=subprocess.DEVNULL,
    )
    return done.returncode


def install(
    repo: Path,
    lab: Path,
    *,
    launchctl: Launchctl,
    notifier: Notifier,
    home: Path | None = None,
    uv: Path | None = None,
    terminal_notifier: Path,
    uid: int | None = None,
    lang: Language = "en",
) -> Path:
    """Write the agent's plist (over whatever is there), load it and send the test
    notification, both in `lang`; raise `ScheduleError` when launchd refuses it or the notifier exits
    non-zero, and then remove the agent again. An agent installed under `OLD_LABEL` is
    removed first. Returns the plist's path."""
    found = uv or (Path(which) if (which := shutil.which("uv")) else None)
    if found is None or not found.is_absolute():
        raise ScheduleError("uv was not found on PATH as an absolute path")
    target = f"gui/{os.getuid() if uid is None else uid}"
    path = agent_path(home if home is not None else Path.home())
    _remove_old(launchctl=launchctl, home=home, uid=uid)
    path.parent.mkdir(parents=True, exist_ok=True)
    (lab / "refresh").mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_bytes(
        plistlib.dumps(
            plist(repo, lab, found, terminal_notifier=terminal_notifier, lang=lang)
        )
    )
    try:
        launchctl(["bootout", f"{target}/{LABEL}"])  # "not loaded" is the usual answer
        # A disabled label cannot be bootstrapped, so it is enabled first.
        if launchctl(["enable", f"{target}/{LABEL}"]) != 0:
            raise ScheduleError(f"launchctl enable {target}/{LABEL} failed")
        if launchctl(["bootstrap", target, str(path)]) != 0:
            raise ScheduleError(f"launchctl bootstrap {target} {path} failed")
        if (status := notifier(TEST_TEXT[lang])) != 0:
            raise ScheduleError(
                f"the test notification failed (terminal-notifier exited {status}): allow"
                " its notifications in System Settings > Notifications, then install again"
            )
    except BaseException:
        uninstall(launchctl=launchctl, home=home, uid=uid)
        raise
    return path


def uninstall(
    *, launchctl: Launchctl, home: Path | None = None, uid: int | None = None
) -> bool:
    """Unload the agent and remove its plist, and one under `OLD_LABEL`; nothing to do
    (False) when neither is there."""
    old = _remove_old(launchctl=launchctl, home=home, uid=uid)
    target = f"gui/{os.getuid() if uid is None else uid}"
    path = agent_path(home if home is not None else Path.home())
    launchctl(["bootout", f"{target}/{LABEL}"])
    if not path.is_file():
        return old
    path.unlink()
    return True


def _remove_old(*, launchctl: Launchctl, home: Path | None, uid: int | None) -> bool:
    """Unload and remove an agent installed under `OLD_LABEL`, when its plist is there."""
    path = agent_path(home if home is not None else Path.home(), OLD_LABEL)
    if not path.is_file():
        return False
    launchctl(["bootout", f"gui/{os.getuid() if uid is None else uid}/{OLD_LABEL}"])
    path.unlink()
    return True


# --- command line ---------------------------------------------------------------------------


def _command(args: argparse.Namespace) -> int:
    if args.action == "uninstall":
        removed = uninstall(launchctl=launchctl)
        print(
            "genokeeper schedule: removed"
            if removed
            else "genokeeper schedule: not installed"
        )
        return 0
    tool = shutil.which("terminal-notifier")
    if tool is None:
        print(
            "genokeeper schedule: terminal-notifier is not installed, and without it no"
            " notification can be proven (brew install terminal-notifier)",
            file=sys.stderr,
        )
        return 1
    found = Path(tool)
    try:
        path = install(
            load_paths().repo,
            load_paths().lab,
            launchctl=launchctl,
            notifier=lambda text: notify(text, terminal_notifier=found),
            terminal_notifier=found,
            lang=language(),
        )
    except ScheduleError as e:
        print(f"genokeeper schedule: {e}", file=sys.stderr)
        return 1
    print(f"genokeeper schedule: installed {path}")
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `schedule` sub-command: `install` and `uninstall`."""
    s = subparsers.add_parser(
        "schedule", help="install or remove the daily update check (a user LaunchAgent)"
    )
    actions = s.add_subparsers(dest="action", required=True)
    actions.add_parser("install", help="write, load and test the agent")
    actions.add_parser("uninstall", help="unload and remove the agent")
    s.set_defaults(handler=_command)
