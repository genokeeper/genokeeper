"""genokeeper.lang: which language user-facing text is written in."""

import pytest

from genokeeper.lang import ENV, language


def test_english_without_a_choice(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(ENV, raising=False)
    assert language() == "en"


def test_the_environment_chooses_when_no_option_does(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(ENV, "ru")
    assert (language(), language("en")) == ("ru", "en")


def test_another_value_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(ENV, "de")
    with pytest.raises(ValueError, match="unknown language 'de': choose one of en, ru"):
        language()
