"""docs/demo/dashboard.html, the README's demo page: its inline script and style still match the
hashes its Content-Security-Policy allows. An edit inside either (a rename did it once) makes the
browser block the script, and the page shows no findings."""

import re
from pathlib import Path

from genokeeper.dashboard import csp_sha256

PAGE = Path(__file__).resolve().parents[1] / "docs" / "demo" / "dashboard.html"


def _one(page: str, tag: str) -> str:
    """The text of the page's one `<tag>` element, which the page writes without attributes."""
    parts = page.split(f"<{tag}>")
    assert len(parts) == 2, f"{len(parts) - 1} <{tag}> elements"
    return parts[1].split(f"</{tag}>", 1)[0]


def test_the_demo_page_matches_its_csp() -> None:
    html = PAGE.read_text(encoding="utf-8")
    csp = re.search(r'http-equiv="Content-Security-Policy" content="([^"]+)"', html)
    script = _one(html, "script")
    style = _one(html, "style")
    assert csp
    assert f"script-src 'sha256-{csp_sha256(script)}'" in csp[1]
    assert f"style-src 'sha256-{csp_sha256(style)}'" in csp[1]
