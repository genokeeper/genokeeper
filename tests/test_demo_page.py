"""docs/demo/dashboard.html, the README's demo page: its inline script and style still match the
hashes its Content-Security-Policy allows. An edit inside either (a rename did it once) makes the
browser block the script, and the page shows no findings."""

import re
from pathlib import Path

from genokeeper.dashboard import _sha256  # pyright: ignore[reportPrivateUsage]

PAGE = Path(__file__).resolve().parents[1] / "docs" / "demo" / "dashboard.html"


def test_the_demo_page_matches_its_csp() -> None:
    html = PAGE.read_text(encoding="utf-8")
    csp = re.search(r'http-equiv="Content-Security-Policy" content="([^"]+)"', html)
    script = re.search(r"<script>(.*?)</script>", html, re.DOTALL)
    style = re.search(r"<style>(.*?)</style>", html, re.DOTALL)
    assert csp and script and style
    assert f"script-src 'sha256-{_sha256(script[1])}'" in csp[1]
    assert f"style-src 'sha256-{_sha256(style[1])}'" in csp[1]
