"""SMN1/SMN2 from the two callers (`genokeeper.smn`): their HG002 outputs (tests/fixtures/hg002)
as the grading reads them, and the two container runs, whose container is replaced by a stub
that checks the arguments and writes the fixture into the output directory."""

import gzip
import json
import shutil
from pathlib import Path

import pytest

from genokeeper import smn
from genokeeper.smn import (
    ParascopyCall,
    SmnCaller,
    parascopy_data,
    parascopy_depth,
    read_parascopy,
    read_smn_caller,
    run_smn_caller,
)

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
OUT_MOUNTS = ("/cram", "/ref", "/out")
SMN_JSON = FIXTURES / "hg002.smncaller.json"
PARALOG = FIXTURES / "hg002.parascopy.paralog.bed.txt"
SMN1_ROW = "chr5\t70947243\t70956598\tHG002\tPASS\t2\t156\t"


def test_read_smn_caller_hg002() -> None:
    assert read_smn_caller(SMN_JSON) == SmnCaller(2, 2, 0, "PASS:Majority")


def test_read_smn_caller_none_values(tmp_path: Path) -> None:
    path = tmp_path / "smn.json"
    path.write_text(
        json.dumps({"sample": {"SMN1": None, "SMN2": "None", "g27134TG_CN": "None"}})
    )
    assert read_smn_caller(path) == SmnCaller(None, None, None, None)


def test_read_parascopy_hg002() -> None:
    assert read_parascopy(PARALOG) == ParascopyCall(2, 2, 156.0, True)


def test_read_parascopy_gzip(tmp_path: Path) -> None:
    path = tmp_path / "res.paralog.bed.gz"
    path.write_bytes(gzip.compress(PARALOG.read_bytes()))
    assert read_parascopy(path) == ParascopyCall(2, 2, 156.0, True)


def test_read_parascopy_unknown(tmp_path: Path) -> None:
    text = PARALOG.read_text()
    assert SMN1_ROW in text
    unknown = text.replace(
        SMN1_ROW, "chr5\t70947243\t70956598\tHG002\tNoReliable\t?\t?\t"
    )
    path = tmp_path / "res.paralog.bed"
    path.write_text(unknown)
    assert read_parascopy(path) == ParascopyCall(None, 2, None, False)
    # SMN2's row is absent: SMN2 is not reported.
    rows = [r for r in text.splitlines() if not r.startswith("chr5\t70071818\t")]
    path.write_text("\n".join(rows) + "\n")
    assert read_parascopy(path) == ParascopyCall(2, None, 156.0, True)


class _Fake:
    """Stands for `run_container`: records each call and writes what SMNCopyNumberCaller and
    Parascopy wrote on HG002."""

    def __init__(self) -> None:
        self.calls: list[
            tuple[str, list[str], list[tuple[Path, str, bool]], str | None]
        ] = []

    def __call__(
        self,
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        self.calls.append((name, list(args), list(mounts), label))
        (out,) = (host for host, target, ro in mounts if target == "/out" and not ro)
        if name == "smncaller":
            shutil.copy(SMN_JSON, out / "smn.json")
        elif args[1] == "cn-using":
            (out / "cn").mkdir()
            (out / "cn" / "res.paralog.bed.gz").write_bytes(
                gzip.compress(PARALOG.read_bytes())
            )


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> _Fake:
    stub = _Fake()
    monkeypatch.setattr(smn, "run_container", stub)
    return stub


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    cram, fasta, out = (
        tmp_path / "in" / "x.cram",
        tmp_path / "ref" / "g.fa",
        tmp_path / "out",
    )
    for path in (cram, fasta):
        path.parent.mkdir()
        path.write_bytes(b"")
    out.mkdir()
    return cram, fasta, out


def _lab(tmp_path: Path) -> Path:
    lab = tmp_path / "lab"
    table = lab / "parascopy" / "homology_table" / "GRCh38.bed.gz"
    model = lab / "parascopy" / "models_GRCh38_1KGP" / "EUR" / "SMN1.gz"
    for path in (table, model):
        path.parent.mkdir(parents=True)
        path.write_bytes(b"")
    return lab


def test_run_smn_caller_args(tmp_path: Path, fake: _Fake) -> None:
    cram, fasta, out = _inputs(tmp_path)
    assert run_smn_caller(cram, fasta, out, threads=3) == out / "smn.json"
    ((name, args, mounts, label),) = fake.calls
    assert (name, label) == ("smncaller", "dna:smn")
    assert args[0] == "smn_caller.py"
    assert args[args.index("--genome") + 1] == "38"
    assert args[args.index("--reference") + 1] == "/ref/g.fa"
    assert args[args.index("--prefix") + 1] == "smn"
    assert args[args.index("--outDir") + 1] == "/out"
    assert args[args.index("--threads") + 1] == "3"
    manifest = args[args.index("--manifest") + 1]
    assert manifest == "/out/manifest.txt"
    assert (out / "manifest.txt").read_text() == "/cram/x.cram\n"
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/out", False),
    }


def test_parascopy_depth_args(tmp_path: Path, fake: _Fake) -> None:
    """The depth step every Parascopy copy-number step reads (3a's SMN1 model and 3d's
    loci), behind `dna:smn`."""
    cram, fasta, out = _inputs(tmp_path)
    assert parascopy_depth(cram, fasta, out, threads=3) == out / "depth" / "depth.csv"
    ((name, args, mounts, label),) = fake.calls
    assert (name, label) == ("parascopy", "dna:smn")
    assert args == [
        "parascopy", "depth", "-i", "/cram/x.cram::sample", "-f", "/ref/g.fa",
        "-@", "3", "-g", "hg38", "-o", "/out/depth",
    ]  # fmt: skip
    assert {(t, ro) for _, t, ro in mounts} == {(m, m != "/out") for m in OUT_MOUNTS}


def test_parascopy_data_in_the_lab(tmp_path: Path) -> None:
    lab = _lab(tmp_path)
    table, models = parascopy_data(lab)
    assert table == lab / "parascopy" / "homology_table" / "GRCh38.bed.gz"
    assert models == lab / "parascopy" / "models_GRCh38_1KGP" / "EUR"


def test_parascopy_data_missing_refuses(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match=r"scripts/bootstrap\.sh"):
        parascopy_data(tmp_path)
    lab = _lab(tmp_path / "half")
    (lab / "parascopy" / "models_GRCh38_1KGP" / "EUR" / "SMN1.gz").unlink()
    with pytest.raises(FileNotFoundError, match=r"SMN1\.gz.*scripts/bootstrap\.sh"):
        parascopy_data(lab)
