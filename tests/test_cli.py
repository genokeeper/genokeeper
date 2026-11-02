import pytest

from genokeeper.cli import main


def test_version_prints_package_version(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--version"]) == 0
    assert capsys.readouterr().out.startswith("genokeeper ")
