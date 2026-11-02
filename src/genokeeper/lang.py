"""The language user-facing text is written in: English by default, Russian on request.

A command with a `--lang` option passes it here; without one, `$GENOKEEPER_LANG` decides, and
without that, English. Every user-facing text that has a Russian version reads its language
from `language()`, so one setting switches the dashboard, the explanations and the messages.
"""

from typing import Literal

from genokeeper.env import getenv

type Language = Literal["en", "ru"]
LANGUAGES: tuple[Language, ...] = ("en", "ru")
ENV = "GENOKEEPER_LANG"  # `$GENOMELAB_LANG` before the rename (`genokeeper.env`)


def language(chosen: str | None = None) -> Language:
    """`chosen` (a command's `--lang`), else `$GENOKEEPER_LANG`, else English.

    Raises `ValueError` for any other value, naming the ones there are."""
    value = chosen or getenv("LANG") or "en"
    if value == "en":
        return "en"
    if value == "ru":
        return "ru"
    raise ValueError(
        f"unknown language {value!r}: choose one of {', '.join(LANGUAGES)}"
    )
