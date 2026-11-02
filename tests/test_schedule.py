"""The LaunchAgent that runs the daily check, and the notifier it calls. `launchctl` and the
notifier are fakes and the home is a tmp directory: nothing here touches the real
`~/Library/LaunchAgents` or shows a notification."""

import plistlib
import subprocess
from pathlib import Path
from typing import Any

import pytest

from genokeeper import cli, schedule
from genokeeper.schedule import LABEL, OLD_LABEL, ScheduleError

UID = 501
UV = Path("/opt/homebrew/bin/uv")
TERMINAL_NOTIFIER = Path("/opt/homebrew/bin/terminal-notifier")
HOSTILE = '1.0" & do shell script "x'


def git(*args: str, cwd: Path) -> None:
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=t",
            "-c",
            "user.email=t@t",
            "-c",
            "commit.gpgsign=false",
            *args,
        ],
        cwd=cwd,
        check=True,
        capture_output=True,
    )


@pytest.fixture
def checkouts(tmp_path: Path) -> tuple[Path, Path]:
    """A repository and a worktree of it: (main checkout, worktree)."""
    main = (tmp_path / "main").resolve()
    main.mkdir()
    git("init", "-q", cwd=main)
    git("commit", "-q", "--allow-empty", "-m", "x", cwd=main)
    worktree = tmp_path / "work" / "tree"
    git("worktree", "add", "-q", str(worktree), "-b", "side", cwd=main)
    return main, worktree.resolve()


class Launchd:
    """A `launchctl` that behaves as the real one does for these three verbs: bootstrapping
    a loaded or a disabled label fails, booting out an unloaded one fails, and enabling
    works whether or not the label is loaded."""

    def __init__(self, *, refuse: str | None = None, disabled: bool = False) -> None:
        self.calls: list[list[str]] = []
        self.loaded = False
        self.refuse = refuse
        self.disabled = disabled

    def __call__(self, args: list[str]) -> int:
        self.calls.append(args)
        verb = args[0]
        if verb == self.refuse:
            return 1
        if verb == "bootout":
            was, self.loaded = self.loaded, False
            return 0 if was else 3
        if verb == "bootstrap":
            if self.loaded or self.disabled:
                return 5
            self.loaded = True
            return 0
        self.disabled = False
        return 0


class Notifier:
    def __init__(self, status: int = 0) -> None:
        self.status = status
        self.texts: list[str] = []

    def __call__(self, text: str) -> int:
        self.texts.append(text)
        return self.status


def install(
    main: Path, tmp_path: Path, launchd: Launchd, notifier: Notifier, **kw: Any
) -> Path:
    return schedule.install(
        main,
        tmp_path / "lab",
        launchctl=launchd,
        notifier=notifier,
        home=tmp_path / "home",
        uv=UV,
        terminal_notifier=TERMINAL_NOTIFIER,
        uid=UID,
        **kw,
    )


def test_plist_has_absolute_paths_and_the_main_checkout(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, worktree = checkouts
    lab = tmp_path / "lab"
    agent = schedule.plist(worktree, lab, UV, terminal_notifier=TERMINAL_NOTIFIER)
    assert agent["Label"] == LABEL == "local.genokeeper.check-updates"
    assert agent["ProgramArguments"] == [
        str(UV),
        "run",
        "--project",
        str(main),
        "--frozen",
        "genokeeper",
        "check-updates",
        "--notify",
    ]
    assert agent["StartCalendarInterval"] == {"Hour": 10, "Minute": 0}
    assert agent["RunAtLoad"] is True
    log = str(lab / "refresh" / "check.log")
    assert (agent["StandardOutPath"], agent["StandardErrorPath"]) == (log, log)
    assert agent["EnvironmentVariables"] == {
        "GENOKEEPER_HOME": str(lab),
        "GENOKEEPER_NOTIFIER": str(TERMINAL_NOTIFIER),
        "GENOKEEPER_LANG": "en",
    }
    # A plist loads from a file: it must survive a round trip, and name nothing relative.
    assert plistlib.loads(plistlib.dumps(agent)) == agent
    absolute = [agent["ProgramArguments"][0], agent["ProgramArguments"][3], log]
    assert all(Path(p).is_absolute() for p in absolute)
    assert str(worktree) not in repr(agent)
    # The main checkout itself is its own project.
    agent = schedule.plist(main, lab, UV, terminal_notifier=TERMINAL_NOTIFIER)
    assert agent["ProgramArguments"][3] == str(main)


def test_plist_refuses_a_directory_outside_git(tmp_path: Path) -> None:
    with pytest.raises(ScheduleError):
        schedule.plist(
            tmp_path, tmp_path / "lab", UV, terminal_notifier=TERMINAL_NOTIFIER
        )


def test_install_writes_loads_and_tests(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    launchd, notifier = Launchd(), Notifier()
    path = install(main, tmp_path, launchd, notifier)
    assert path == tmp_path / "home" / "Library" / "LaunchAgents" / f"{LABEL}.plist"
    assert plistlib.loads(path.read_bytes()) == schedule.plist(
        main, tmp_path / "lab", UV, terminal_notifier=TERMINAL_NOTIFIER
    )
    assert launchd.calls == [
        ["bootout", f"gui/{UID}/{LABEL}"],
        ["enable", f"gui/{UID}/{LABEL}"],
        ["bootstrap", f"gui/{UID}", str(path)],
    ]
    assert len(notifier.texts) == 1  # the test notification, after the agent is loaded
    assert notifier.texts == [schedule.TEST_TEXT["en"]]


def test_install_in_russian_tests_and_runs_the_agent_in_russian(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    notifier = Notifier()
    path = install(main, tmp_path, Launchd(), notifier, lang="ru")
    assert notifier.texts == ["Проверка обновлений включена: это тестовое уведомление."]
    agent = plistlib.loads(path.read_bytes())
    assert agent["EnvironmentVariables"]["GENOKEEPER_LANG"] == "ru"


def test_install_twice_boots_out_first(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    launchd, notifier = Launchd(), Notifier()
    install(main, tmp_path, launchd, notifier)
    install(main, tmp_path, launchd, notifier)  # bootstrap would fail on a loaded label
    verbs = [call[0] for call in launchd.calls]
    assert verbs == ["bootout", "enable", "bootstrap"] * 2
    assert launchd.loaded and len(notifier.texts) == 2


def test_install_fails_when_the_test_notification_fails(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    """A failed install leaves no agent behind: it would run daily into the same failure."""
    main, _ = checkouts
    launchd = Launchd()
    with pytest.raises(ScheduleError, match=r"exited 1\): allow its notifications"):
        install(main, tmp_path, launchd, Notifier(status=1))
    assert not launchd.loaded
    assert launchd.calls[-1] == ["bootout", f"gui/{UID}/{LABEL}"]
    assert not schedule.agent_path(tmp_path / "home").exists()


@pytest.mark.parametrize("verb", ["bootstrap", "enable"])
def test_install_fails_when_launchd_refuses(
    checkouts: tuple[Path, Path], tmp_path: Path, verb: str
) -> None:
    main, _ = checkouts
    notifier, launchd = Notifier(), Launchd(refuse=verb)
    with pytest.raises(ScheduleError, match=verb):
        install(main, tmp_path, launchd, notifier)
    assert notifier.texts == []  # no test notification for an agent that did not load
    # A refusal after the plist is written removes the agent again.
    assert not launchd.loaded
    assert launchd.calls[-1] == ["bootout", f"gui/{UID}/{LABEL}"]
    assert not schedule.agent_path(tmp_path / "home").exists()


def test_install_re_enables_a_disabled_label(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    """A disabled label cannot be bootstrapped until it is enabled: enable comes first."""
    main, _ = checkouts
    launchd = Launchd(disabled=True)
    install(main, tmp_path, launchd, Notifier())
    assert launchd.loaded and not launchd.disabled


def test_install_needs_an_absolute_uv(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    with pytest.raises(ScheduleError, match="uv"):
        schedule.install(
            main,
            tmp_path / "lab",
            launchctl=Launchd(),
            notifier=Notifier(),
            home=tmp_path / "home",
            uv=Path("uv"),
            terminal_notifier=TERMINAL_NOTIFIER,
            uid=UID,
        )


def test_uninstall_removes_the_agent(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    launchd = Launchd()
    path = install(main, tmp_path, launchd, Notifier())
    assert schedule.uninstall(launchctl=launchd, home=tmp_path / "home", uid=UID)
    assert not path.exists() and not launchd.loaded
    assert launchd.calls[-1] == ["bootout", f"gui/{UID}/{LABEL}"]


def test_uninstall_when_absent_is_quiet(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    launchd = Launchd()
    home = tmp_path / "home"
    assert schedule.uninstall(launchctl=launchd, home=home, uid=UID) is False
    assert not home.exists()  # nothing was created to remove it
    assert launchd.calls == [["bootout", f"gui/{UID}/{LABEL}"]]
    out, err = capsys.readouterr()
    assert (out, err) == ("", "")


def old_agent(tmp_path: Path) -> Path:
    """An agent installed before the rename from genomelab, under its label then."""
    path = schedule.agent_path(tmp_path / "home", OLD_LABEL)
    path.parent.mkdir(parents=True)
    path.write_bytes(plistlib.dumps({"Label": OLD_LABEL}))
    return path


def test_install_removes_the_agent_installed_before_the_rename(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    old = old_agent(tmp_path)
    launchd = Launchd()
    path = install(main, tmp_path, launchd, Notifier())
    assert not old.exists() and path.is_file()
    assert launchd.calls[0] == ["bootout", f"gui/{UID}/{OLD_LABEL}"]


def test_uninstall_removes_the_agent_installed_before_the_rename(
    tmp_path: Path,
) -> None:
    old = old_agent(tmp_path)
    launchd = Launchd()
    assert schedule.uninstall(launchctl=launchd, home=tmp_path / "home", uid=UID)
    assert not old.exists()
    assert ["bootout", f"gui/{UID}/{OLD_LABEL}"] in launchd.calls


def test_hand_edited_plist_is_rewritten_on_install(
    checkouts: tuple[Path, Path], tmp_path: Path
) -> None:
    main, _ = checkouts
    launchd = Launchd()
    path = install(main, tmp_path, launchd, Notifier())
    edited = plistlib.loads(path.read_bytes())
    edited["ProgramArguments"] = ["/bin/sh", "-c", "echo hand-edited"]
    edited["StartInterval"] = 5
    path.write_bytes(plistlib.dumps(edited))
    install(main, tmp_path, launchd, Notifier())
    assert plistlib.loads(path.read_bytes()) == schedule.plist(
        main, tmp_path / "lab", UV, terminal_notifier=TERMINAL_NOTIFIER
    )
    path.write_text("not a plist at all")
    install(main, tmp_path, launchd, Notifier())
    assert plistlib.loads(path.read_bytes())["Label"] == LABEL


class Run:
    """A stand-in for `subprocess.run` that records the argument vector."""

    def __init__(self) -> None:
        self.argv: list[str] = []

    def __call__(
        self, argv: list[str], **kwargs: Any
    ) -> subprocess.CompletedProcess[bytes]:
        self.argv = argv
        assert kwargs["check"] is False
        return subprocess.CompletedProcess(argv, 7)


def test_notify_passes_text_as_argv() -> None:
    """A `"` in the text never reaches a script: terminal-notifier gets it as the value of
    `-message`."""
    run = Run()
    status = schedule.notify(HOSTILE, terminal_notifier=TERMINAL_NOTIFIER, run=run)
    assert status == 7  # the notifier's own exit status
    assert run.argv == [
        str(TERMINAL_NOTIFIER),
        "-title",
        schedule.TITLE,
        "-message",
        HOSTILE,
    ]

    # Without terminal-notifier there is no notifier: osascript can exit 0 and show nothing.
    run.argv = []
    assert schedule.notify(HOSTILE, environ={}, run=run) != 0
    assert run.argv == []


def test_notify_shows_the_notifiers_own_reason(
    capsys: pytest.CaptureFixture[str],
) -> None:
    reason = (
        b"Could not request notification permission: Notifications are not allowed\n"
    )

    def refused(argv: list[str], **_: Any) -> subprocess.CompletedProcess[bytes]:
        return subprocess.CompletedProcess(argv, 3, b"", reason)

    assert schedule.notify("t", terminal_notifier=TERMINAL_NOTIFIER, run=refused) == 3
    assert capsys.readouterr().err == reason.decode()


def test_notify_finds_terminal_notifier_in_the_agents_environment() -> None:
    run = Run()
    schedule.notify(
        "t", environ={schedule.NOTIFIER_ENV: str(TERMINAL_NOTIFIER)}, run=run
    )
    assert run.argv[0] == str(TERMINAL_NOTIFIER)


def test_notify_finds_terminal_notifier_under_its_name_before_the_rename() -> None:
    """An agent installed before the rename sets `GENOMELAB_NOTIFIER`."""
    run = Run()
    schedule.notify(
        "t", environ={"GENOMELAB_NOTIFIER": str(TERMINAL_NOTIFIER)}, run=run
    )
    assert run.argv[0] == str(TERMINAL_NOTIFIER)


def test_the_commands_are_registered(
    checkouts: tuple[Path, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    main, _ = checkouts
    launchd, notifier = Launchd(), Notifier()
    monkeypatch.setattr(schedule, "launchctl", launchd)
    monkeypatch.setattr(schedule, "notify", lambda text, **kw: notifier(text))
    monkeypatch.setattr(schedule.Path, "home", lambda: tmp_path / "home")
    monkeypatch.setattr(
        schedule.shutil, "which", lambda name: f"/opt/homebrew/bin/{name}"
    )
    paths = type("P", (), {"repo": main, "lab": tmp_path / "lab"})()
    monkeypatch.setattr(schedule, "load_paths", lambda: paths)
    monkeypatch.delenv("GENOKEEPER_LANG", raising=False)
    assert cli.main(["schedule", "install"]) == 0
    assert (tmp_path / "home" / "Library" / "LaunchAgents" / f"{LABEL}.plist").is_file()
    monkeypatch.setenv("GENOKEEPER_LANG", "ru")
    assert cli.main(["schedule", "install"]) == 0
    assert notifier.texts == [schedule.TEST_TEXT["en"], schedule.TEST_TEXT["ru"]]
    assert cli.main(["schedule", "uninstall"]) == 0
    assert cli.main(["schedule", "uninstall"]) == 0
    assert "not installed" in capsys.readouterr().out
    # Without terminal-notifier nothing is installed: no notification could be proven.
    monkeypatch.setattr(
        schedule.shutil,
        "which",
        lambda name: None if "notifier" in name else f"/opt/homebrew/bin/{name}",
    )
    calls = len(launchd.calls)
    assert cli.main(["schedule", "install"]) == 1
    assert "brew install terminal-notifier" in capsys.readouterr().err
    assert len(launchd.calls) == calls
    assert not (
        tmp_path / "home" / "Library" / "LaunchAgents" / f"{LABEL}.plist"
    ).exists()
