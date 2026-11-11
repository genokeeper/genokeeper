"""A packaged HG002 demo (`demo.package.package`'s `out/`) on a small two-run store, and a
local HTTP server that serves it as the demo's release does."""

import http.server
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from fixtures.store_mini import Row, seed
from genokeeper import refresh_log
from genokeeper.demo import STORE_FILE, demo_paths
from genokeeper.demo.package import package
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus
from genokeeper.store.db import connect

RUN1, RUN2 = "20261012T120000Z-aaaaaa", "20261109T120000Z-bbbbbb"
GENE = "GENE1"


def seeded(lab: Path, before: list[Row], after: list[Row]) -> Path:
    """The demo root under `lab`: its store (returned) with RUN1's `before` and RUN2's
    `after`, and the refresh-log record that made RUN2."""
    root = demo_paths(lab)
    store = root.reports / STORE_FILE
    with connect(store) as con:
        seed(con, before, run_id=RUN1)
        seed(con, after, run_id=RUN2)
    at = datetime(2026, 11, 10, 9, tzinfo=UTC)
    clinvar = SourceOutcome("clinvar", "20261011", "20261108", SourceStatus.NEW)
    record = RefreshRecord(at, at, "complete", (clinvar,), RUN1, None, RUN2)
    refresh_log.append(root.repo / "lab", record)
    return store


def packaged(tmp: Path) -> Path:
    """`package`'s `out/` under `tmp`: two runs, one finding whose evidence moved between
    them, and the refresh-log record that made the second."""
    home = tmp / "home"
    lab = home / "genome-lab"
    rows = [Row("clinvar:2", gene="GENE2", pos=2000, validity="moderate")]
    seeded(
        lab,
        [Row("clinvar:1", gene=GENE, validity="limited"), *rows],
        [Row("clinvar:1", gene=GENE, validity="strong"), *rows],
    )
    return package(lab, tmp / "out", home=home, private_ids=None)


@contextmanager
def serve(folder: Path, *, cut_store_at: int | None = None) -> Iterator[str]:
    """`folder` over HTTP on 127.0.0.1, as a base URL. With `cut_store_at`, the store's
    response promises its whole length and the connection closes after that many bytes."""

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, directory=str(folder), **kwargs)

        def log_message(self, format: str, *args: Any) -> None:
            pass

        def do_GET(self) -> None:
            if cut_store_at is None or self.path != f"/{STORE_FILE}":
                return super().do_GET()
            data = (folder / STORE_FILE).read_bytes()
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data[:cut_store_at])
            self.close_connection = True

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_address[1]}/"
    finally:
        server.shutdown()
        server.server_close()
