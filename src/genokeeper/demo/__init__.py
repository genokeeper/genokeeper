"""The HG002 demo: GIAB's HG002 run through genokeeper, packaged for the public site and for
`genokeeper demo`.

The maintainer builds it under the lab, in `DEMO_ROOT` (`demo_paths`: its own `reports/`
and `lab/`, never the repository's or the real lab's), and `python -m genokeeper.demo package
OUT` writes what is published: the store (`STORE_FILE`), the refresh-log record
(`REFRESH_LOG_FILE`), the site's pages (`SITE_FILES`) and `MANIFEST`, which names the
version, the store's schema version, the runs and each file's sha256.
"""

from pathlib import Path

from genokeeper.config import Paths

DEMO_ROOT = Path("demo/hg002")
MANIFEST = "manifest.json"
STORE_FILE = "genomelab.duckdb"
REFRESH_LOG_FILE = "refresh-log.jsonl"
SITE_FILES = ("index.html", "dashboard.html", "calibration.html")


def demo_paths(lab: Path) -> Paths:
    """The demo's paths: its root `lab/DEMO_ROOT` with its own `raw/` and `reports/`, beside
    the real lab whose reference data and tools it runs with."""
    repo = lab / DEMO_ROOT
    return Paths(repo=repo, raw=repo / "raw", reports=repo / "reports", lab=lab)
