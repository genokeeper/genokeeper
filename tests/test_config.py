import dataclasses
from pathlib import Path

import pytest

from genokeeper import config
from genokeeper.config import (
    DEMO_MARKER,
    AfGuardThresholds,
    CnvThresholds,
    DeletionThresholds,
    DosageThresholds,
    Grading,
    Hhv6Thresholds,
    InsertionThresholds,
    LpaThresholds,
    ParalogThresholds,
    QcThresholds,
    RohThresholds,
    SexThresholds,
    load_grading,
    load_paths,
)
from genokeeper.repeats.loci import Rule
from genokeeper.store.db import Validity


def test_grading_defaults() -> None:
    grading = load_grading()
    assert grading.p_high == 1e-3
    assert grading.p_low == 1e-2
    # GRCh38 pseudoautosomal regions, 1-based inclusive (GIAB v3.6 XY/GRCh38_chr{X,Y}_PAR).
    assert grading.par == {
        "chrX": [(10001, 2781479), (155701383, 156030895)],
        "chrY": [(10001, 2781479), (56887903, 57217415)],
    }


def test_the_dataclass_defaults_are_the_shipped_values() -> None:
    shipped = load_grading()
    built_in = Grading(
        ba1=shipped.ba1, bs1=shipped.bs1, paralogs=shipped.paralogs, lpa=shipped.lpa
    )
    tables = {"par", "bs1_by_gene", "strata"}  # data tables loaded from the file
    for field in dataclasses.fields(Grading):
        if field.name in tables:
            continue
        mine, theirs = getattr(built_in, field.name), getattr(shipped, field.name)
        if field.name == "repeats":  # its `loci` table is data too
            theirs = dataclasses.replace(theirs, loci={})
        assert mine == theirs, field.name


SHIPPED = Path(__file__).parents[1] / "config" / "grading.toml"


def _variant(tmp_path: Path, old: str, new: str) -> Path:
    """The shipped config with one line changed."""
    text = SHIPPED.read_text()
    assert text.count(old) == 1, old
    config = tmp_path / "grading.toml"
    config.write_text(text.replace(old, new))
    return config


def test_grading_rejects_unknown_keys(tmp_path: Path) -> None:
    config = tmp_path / "grading.toml"
    config.write_text(SHIPPED.read_text() + "\n[extra]\np_hihg = 1e-4\n")
    with pytest.raises(ValueError, match="extra"):
        load_grading(config)


def test_grading_rejects_a_misspelt_key_in_a_known_section(tmp_path: Path) -> None:
    config = _variant(tmp_path, "\nbs1 = ", "\nbs_1 = ")
    with pytest.raises(ValueError, match=r"\[frequency\].*bs_1"):
        load_grading(config)


@pytest.mark.parametrize(
    ("old", "new", "key"),
    [
        ("min_dp = 10", 'min_dp = "10"', "min_dp"),
        ("min_dp = 10", "min_dp = 10.5", "min_dp"),
        ("min_dp = 10", "min_dp = true", "min_dp"),
        ("p_high = 1e-3", 'p_high = "1e-3"', "p_high"),
        ("ba1 = 0.05", "ba1 = true", "ba1"),
        ("[10001, 2781479], [155701383", '[10001, "2781479"], [155701383', "chrX"),
        ("[10001, 2781479], [155701383", "[10001], [155701383", "chrX"),
        ('mhc = "GRCh38_MHC.bed.gz"', "mhc = 5", "mhc"),
    ],
)
def test_grading_rejects_values_of_the_wrong_type(
    tmp_path: Path, old: str, new: str, key: str
) -> None:
    with pytest.raises(TypeError, match=key):
        load_grading(_variant(tmp_path, old, new))


def test_bs1_for_takes_the_gene_value_else_the_fallback() -> None:
    grading = load_grading()
    # Johnston & Biesecker 2025: Table 1 (VCEP), Table 2 (one entity), Table 3 (several).
    assert grading.bs1_for("BRCA1") == 0.0001
    assert grading.bs1_for("TTN") == 0.001299
    assert grading.bs1_for("SCN5A") == 0.00103  # the highest Table 3 value for SCN5A
    assert grading.bs1_for("F5") == grading.bs1 == 0.01
    assert len(grading.bs1_by_gene) == 83


def test_bs1_for_a_recessive_gene_outside_the_table_is_none() -> None:
    grading = load_grading()
    assert grading.bs1_for("CFTR", autosomal_recessive=True) is None
    assert grading.bs1_for("CFTR") == grading.bs1


def test_bs1_for_a_listed_gene_is_its_value_either_way() -> None:
    grading = load_grading()
    assert grading.bs1_for("GAA", autosomal_recessive=True) == 0.005
    assert grading.bs1_for("GAA", autosomal_recessive=False) == 0.005


def test_bs1_by_gene_values_are_type_and_range_checked(tmp_path: Path) -> None:
    with pytest.raises(TypeError, match="TTN"):
        load_grading(_variant(tmp_path, "TTN = 0.001299", 'TTN = "0.001299"'))
    with pytest.raises(ValueError, match="TTN"):
        load_grading(_variant(tmp_path, "TTN = 0.001299", "TTN = 0.2"))


def test_grading_needs_the_par_section(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="missing sections.*par"):
        load_grading(_variant(tmp_path, "[par]", "[unused_par]"))


def test_grading_takes_an_integer_where_a_float_is_due(tmp_path: Path) -> None:
    grading = load_grading(
        _variant(tmp_path, "max_dp_factor = 2.0", "max_dp_factor = 2")
    )
    assert grading.max_dp_factor == 2.0
    assert isinstance(grading.max_dp_factor, float)


def test_lab_defaults_to_home_and_env_overrides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("GENOKEEPER_HOME", raising=False)
    paths = load_paths(tmp_path)
    assert paths.repo == tmp_path
    assert paths.raw == tmp_path / "data" / "raw"
    assert paths.reports == tmp_path / "reports"
    assert paths.lab == Path.home() / "genome-lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path / "lab"))
    assert load_paths(tmp_path).lab == tmp_path / "lab"


def test_grch38_sections_load() -> None:
    grading = load_grading()
    assert grading.strata["segdup"] == "GRCh38_segdups.bed.gz"
    assert all(name.startswith("GRCh38_") for name in grading.strata.values())
    assert grading.qc == QcThresholds(
        min_mean_depth=15.0,
        depth_mapq=20,
        min_fraction_20x=0.80,
        min_mapped_fraction=0.98,
        max_duplicate_fraction=0.10,
        min_median_insert=250.0,
        max_freemix=0.03,
        min_titv=1.9,
        max_titv=2.3,
        max_het_hom=3.0,
    )
    assert grading.sex == SexThresholds(mapq=30)
    # The dataclass defaults are the shipped values, for a Grading built in code.
    built_in = Grading(ba1=0.05, bs1=0.01, paralogs=grading.paralogs, lpa=grading.lpa)
    assert built_in.qc == grading.qc
    assert built_in.sex == grading.sex


@pytest.mark.parametrize(
    ("old", "new", "error", "match"),
    [
        ("max_freemix = 0.03", "", ValueError, r"\[qc\].*max_freemix"),
        ("max_freemix = 0.03", "max_freemix = 1.5", ValueError, "max_freemix"),
        ("depth_mapq = 20", "depth_mapq = 20.5", TypeError, "depth_mapq"),
        ("min_titv = 1.9", "min_titv = 2.4", ValueError, "titv"),
        ("[sex]", "[sexx]", ValueError, r"missing sections.*'sex'"),
    ],
)
def test_qc_and_sex_values_checked(
    tmp_path: Path, old: str, new: str, error: type[Exception], match: str
) -> None:
    with pytest.raises(error, match=match):
        load_grading(_variant(tmp_path, old, new))


def test_repeat_loci_loaded() -> None:
    grading = load_grading()
    repeats = grading.repeats
    eleven = [
        "HTT", "AR", "DMPK", "C9orf72", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A",
        "ATN1", "FXN",
    ]  # fmt: skip
    thirteen = [
        "TCF4", "PABPN1", "FMR1", "AFF2", "CNBP", "ATXN10", "NOP56", "PPP2R2B", "JPH3",
        "ZFHX3", "THAP11", "TBP", "FGF14",
    ]  # fmt: skip
    assert list(repeats.loci) == eleven + thirteen
    htt = repeats.loci["HTT"]
    # STRchive's pathogenic minimum is 36; this design keeps 36-39 as reduced penetrance.
    assert (htt.reduced_penetrance, htt.pathogenic_min) == ((36, 39), 40)
    assert repeats.loci["AR"].chrom == "chrX"
    assert repeats.loci["FXN"].topic == "carrier"
    assert all(repeats.loci[g].evidence == "Definitive" for g in eleven)
    rules = {g: repeats.loci[g].rule for g in ("FMR1", "FGF14", "CNBP", "HTT")}
    assert rules == {
        "FMR1": Rule.X_LINKED,
        "FGF14": Rule.DISTINCT,
        "CNBP": Rule.EXCLUSION,
        "HTT": Rule.STANDARD,
    }
    assert repeats.loci["TCF4"].validity is Validity.LIMITED
    assert repeats.loci["HTT"].validity is None
    assert repeats.loci["PABPN1"].recessive_at == 11
    # ATXN10's intermediate stops below its own pathogenic minimum (STRchive: 33-850).
    assert repeats.loci["ATXN10"].intermediate == (33, 799)
    assert repeats.loci["THAP11"].intermediate is None
    assert grading.exclusion.min_spanning == 3


def test_x_linked_on_an_autosome_refused(tmp_path: Path) -> None:
    config = _variant(
        tmp_path, 'motif = "GGCCCC"', 'motif = "GGCCCC", rule = "x_linked"'
    )
    with pytest.raises(ValueError, match=r"\[repeat_loci\] C9orf72.*chrX"):
        load_grading(config)


def test_exclusion_on_chrx_refused(tmp_path: Path) -> None:
    config = _variant(
        tmp_path,
        'variant_id = "AR", chrom = "chrX", motif = "GCA"',
        'variant_id = "AR", chrom = "chrX", motif = "GCA", rule = "exclusion"',
    )
    with pytest.raises(ValueError, match=r"\[repeat_loci\] AR.*chrX"):
        load_grading(config)
    distinct = _variant(
        tmp_path,
        'variant_id = "AR", chrom = "chrX", motif = "GCA"',
        'variant_id = "AR", chrom = "chrX", motif = "GCA", rule = "distinct"',
    )
    with pytest.raises(ValueError, match=r"\[repeat_loci\] AR.*chrX"):
        load_grading(distinct)


def test_unknown_rule_refused(tmp_path: Path) -> None:
    config = _variant(tmp_path, 'motif = "GGCCCC"', 'motif = "GGCCCC", rule = "loose"')
    with pytest.raises(ValueError, match=r"\[repeat_loci\] C9orf72.*loose"):
        load_grading(config)


def test_min_spanning_floor(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"\[exclusion\] min_spanning"):
        load_grading(_variant(tmp_path, "min_spanning = 3", "min_spanning = 2"))
    with pytest.raises(ValueError, match=r"missing sections.*'exclusion'"):
        load_grading(_variant(tmp_path, "[exclusion]", "[unused_exclusion]"))


def test_recessive_at_outside_the_gap_refused(tmp_path: Path) -> None:
    for value in (10, 12):
        config = _variant(tmp_path, "recessive_at = 11", f"recessive_at = {value}")
        with pytest.raises(ValueError, match=r"\[repeat_loci\] PABPN1.*recessive_at"):
            load_grading(config)


def test_repeat_locus_unknown_key_refused(tmp_path: Path) -> None:
    config = _variant(tmp_path, 'motif = "GGCCCC"', 'motif = "GGCCCC", motiff = "x"')
    with pytest.raises(ValueError, match=r"\[repeat_loci\] C9orf72.*motiff"):
        load_grading(config)


def test_deletion_thresholds_loaded() -> None:
    d = load_grading().deletions
    assert d == DeletionThresholds()
    assert (d.min_length, d.long_from, d.exon_margin) == (50, 1000, 20)
    assert (d.max_gnomad_af, d.min_gnomad_an, d.min_reciprocal_overlap) == (
        0.01,
        2000,
        0.5,
    )
    assert (d.max_depth_ratio, d.likely_below) == (0.7, 0.1)


def test_deletion_unknown_key_refused(tmp_path: Path) -> None:
    config = _variant(tmp_path, "long_from = 1000", "long_from = 1000\nlong_form = 2")
    with pytest.raises(ValueError, match=r"\[deletions\].*long_form"):
        load_grading(config)


def test_exon_margin_must_leave_no_intron_a_deletion_fits_in(tmp_path: Path) -> None:
    """Exons at most twice the margin apart merge into one target, so a margin of half
    `min_length` or more lets a kept deletion lie wholly in an intron inside a target, where
    the screen's flank test (`screen.flank_only`) reads it as an exon hit."""
    assert load_grading(_variant(tmp_path, "exon_margin = 20", "exon_margin = 24"))
    with pytest.raises(
        ValueError, match=r"\[deletions\] need 2 x exon_margin < min_length"
    ):
        load_grading(_variant(tmp_path, "exon_margin = 20", "exon_margin = 25"))


def test_insertion_thresholds_loaded(tmp_path: Path) -> None:
    """[insertions] holds only its distance: every other value is [deletions]'s."""
    g = load_grading()
    assert g.insertions == InsertionThresholds(
        **dataclasses.asdict(g.deletions), max_distance=500
    )
    moved = load_grading(_variant(tmp_path, "long_from = 1000", "long_from = 2000"))
    assert (moved.insertions.long_from, moved.insertions.max_distance) == (2000, 500)
    with pytest.raises(ValueError, match=r"\[insertions\].*min_length"):
        load_grading(
            _variant(
                tmp_path, "max_distance = 500", "max_distance = 500\nmin_length = 9"
            )
        )
    with pytest.raises(ValueError, match=r"\[insertions\] max_distance"):
        load_grading(_variant(tmp_path, "max_distance = 500", "max_distance = -1"))
    with pytest.raises(ValueError, match=r"\[insertions\] lacks \['max_distance'\]"):
        load_grading(_variant(tmp_path, "max_distance = 500", ""))


def test_cnv_thresholds_loaded() -> None:
    c = load_grading().cnv
    assert c == CnvThresholds()
    assert (c.min_length, c.band_from, c.large_from, c.bin_size) == (
        1000,
        5000,
        500_000,
        500,
    )


def test_cnv_unknown_key_refused(tmp_path: Path) -> None:
    config = _variant(tmp_path, "bin_size = 500", "bin_size = 500\nbin_sise = 1000")
    with pytest.raises(ValueError, match=r"\[cnv\].*bin_sise"):
        load_grading(config)
    with pytest.raises(ValueError, match=r"\[cnv\] need 0 < min_length < band_from"):
        load_grading(_variant(tmp_path, "band_from = 5000", "band_from = 600000"))
    with pytest.raises(ValueError, match=r"\[cnv\] lacks \['large_from'\]"):
        load_grading(_variant(tmp_path, "large_from = 500000", ""))


def test_paralog_thresholds_loaded(tmp_path: Path) -> None:
    p = load_grading().paralogs
    # `gba_min_reads` has no default: Gauchian's own floor, read with its source.
    assert p == ParalogThresholds(gba_min_reads=5)
    assert (p.min_qual, p.unique_loss_ratio, p.gba_min_reads) == (20.0, (0.3, 0.7), 5)
    with pytest.raises(ValueError, match=r"\[paralogs\] lacks \['gba_min_reads'\]"):
        load_grading(_variant(tmp_path, "gba_min_reads = 5", ""))
    with pytest.raises(ValueError, match=r"\[paralogs\] gba_min_reads"):
        load_grading(_variant(tmp_path, "gba_min_reads = 5", "gba_min_reads = 0"))
    with pytest.raises(ValueError, match=r"\[paralogs\] unique_loss_ratio"):
        load_grading(
            _variant(
                tmp_path,
                "unique_loss_ratio = [0.3, 0.7]",
                "unique_loss_ratio = [0.7, 0.3]",
            )
        )
    with pytest.raises(TypeError, match=r"\[paralogs\] unique_loss_ratio"):
        load_grading(
            _variant(
                tmp_path, "unique_loss_ratio = [0.3, 0.7]", "unique_loss_ratio = 0.5"
            )
        )
    with pytest.raises(ValueError, match=r"\[paralogs\].*min_quall"):
        load_grading(
            _variant(tmp_path, "gba_min_reads = 5", "gba_min_reads = 5\nmin_quall = 1")
        )


def test_lpa_thresholds_loaded(tmp_path: Path) -> None:
    lpa = load_grading().lpa
    # A fraction of HG002's assembly-derived KIV-2 count, not a number of copies.
    assert lpa == LpaThresholds(max_kiv2_error=0.20)
    for bad in ("0", "1.5"):
        with pytest.raises(ValueError, match=r"\[lpa\] max_kiv2_error"):
            load_grading(
                _variant(tmp_path, "max_kiv2_error = 0.20", f"max_kiv2_error = {bad}")
            )
    with pytest.raises(ValueError, match=r"missing sections.*'lpa'"):
        load_grading(_variant(tmp_path, "[lpa]", "[unused_lpa]"))


@pytest.mark.parametrize(
    ("old", "new", "match"),
    [
        ("[strata]", "[unused_strata]", r"missing sections.*'strata'"),
        ("[reliability]", "[unused_reliability]", r"missing sections.*'reliability'"),
        ("p_low = 1e-2", "", r"\[reliability\] lacks \['p_low'\]"),
        ("ba1 = 0.05", "", r"\[frequency\] lacks \['ba1'\]"),
        ("homopolymer = 6", "homopolymer = 0", r"\[predicted\] homopolymer"),
        ("homopolymer = 6", "homopolymer = 1", r"\[predicted\] homopolymer"),
    ],
)
def test_grading_inputs_never_fall_back_silently(
    tmp_path: Path, old: str, new: str, match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        load_grading(_variant(tmp_path, old, new))


def test_identity_section_required(tmp_path: Path) -> None:
    identity = load_grading().identity
    assert (identity.min_concordance, identity.min_compared) == (0.98, 1000)
    with pytest.raises(ValueError, match=r"missing sections.*'identity'"):
        load_grading(_variant(tmp_path, "[identity]", "[unused_identity]"))
    with pytest.raises(ValueError, match=r"\[identity\] min_concordance"):
        load_grading(
            _variant(tmp_path, "min_concordance = 0.98", "min_concordance = 1.2")
        )
    with pytest.raises(ValueError, match=r"\[identity\] min_compared"):
        load_grading(_variant(tmp_path, "min_compared = 1000", "min_compared = 0"))


def test_coverage_section_required(tmp_path: Path) -> None:
    coverage = load_grading().coverage
    assert (coverage.mapq, coverage.min_depth_per_copy, coverage.cds_margin) == (
        20,
        10,
        2,
    )
    assert coverage.min_callable == 0.95
    with pytest.raises(ValueError, match=r"missing sections.*'coverage'"):
        load_grading(_variant(tmp_path, "[coverage]", "[unused_coverage]"))
    with pytest.raises(ValueError, match=r"\[coverage\] min_callable"):
        load_grading(_variant(tmp_path, "min_callable = 0.95", "min_callable = 1.5"))


def test_dosage_section_checked(tmp_path: Path) -> None:
    d = load_grading().dosage
    assert d == DosageThresholds()
    # The design's values, which the calibration measures under and never moves.
    assert (d.tolerance, d.het_low, d.het_high) == (0.15, 0.2, 0.8)
    assert d.min_x_sites == 1000
    assert d.calibrated_on  # the report's date and controls, the statement's version
    for old, new in [
        ("tolerance = 0.15", "tolerance = 0.25"),
        ("het_low = 0.2", "het_low = 0.8"),
        (f"y_absent_max = {d.y_absent_max}", "y_absent_max = 0.35"),
    ]:
        with pytest.raises(ValueError, match=r"\[dosage\]"):
            load_grading(_variant(tmp_path, old, new))
    with pytest.raises(TypeError, match="calibrated_on"):
        load_grading(
            _variant(
                tmp_path,
                f'calibrated_on = "{d.calibrated_on}"',
                "calibrated_on = 1",
            )
        )
    with pytest.raises(ValueError, match=r"unknown keys in \[sex\]"):
        load_grading(
            _variant(tmp_path, "mapq = 30", "mapq = 30\nmale_min_ratio = 0.35")
        )


def test_af_guard_section_required(tmp_path: Path) -> None:
    g = load_grading().af_guard
    assert g == AfGuardThresholds()
    assert g.calibrated_on
    assert g.clonal_genes == (
        "TP53", "CHEK2", "DNMT3A", "TET2", "ASXL1", "PPM1D", "JAK2", "SF3B1", "SRSF2",
    )  # fmt: skip
    with pytest.raises(ValueError, match=r"missing sections.*'af_guard'"):
        load_grading(_variant(tmp_path, "[af_guard]", "[unused_af_guard]"))
    with pytest.raises(ValueError, match=r"\[af_guard\] min_het_vaf"):
        load_grading(
            _variant(tmp_path, f"min_het_vaf = {g.min_het_vaf}", "min_het_vaf = 0.6")
        )
    with pytest.raises(TypeError, match=r"\[af_guard\] clonal_genes"):
        load_grading(_variant(tmp_path, '"TP53", ', "53, "))
    with pytest.raises(ValueError, match=r"\[af_guard\] clonal_genes"):
        load_grading(_variant(tmp_path, '"TP53", ', '"CHEK2", '))


def test_config_rejects_floor_outside_range(tmp_path: Path) -> None:
    """[mito] heteroplasmy_floor: at least 0.02 (mutserve's noise lies below) and below
    min_homoplasmy, else a heteroplasmy could never be graded."""
    line = "\nheteroplasmy_floor = 0.10\n"
    for bad in ("0.019", "0.9"):
        config = _variant(tmp_path, line, f"\nheteroplasmy_floor = {bad}\n")
        with pytest.raises(ValueError, match="heteroplasmy_floor"):
            load_grading(config)
    config = _variant(tmp_path, line, "\nheteroplasmy_floor = 0.02\n")
    assert load_grading(config).mito.heteroplasmy_floor == 0.02


def test_config_hhv6_roh_tables(tmp_path: Path) -> None:
    """[hhv6] (the band around one integrated copy and the background factor) and [roh]
    (the long-run and isodisomy lengths), checked."""
    grading = load_grading()
    assert grading.hhv6 == Hhv6Thresholds()
    assert (grading.hhv6.min_ratio, grading.hhv6.max_ratio) == (0.5, 1.5)
    assert grading.hhv6.background_factor == 3.0
    assert grading.roh == RohThresholds()
    assert (grading.roh.min_run, grading.roh.isodisomy_min) == (5_000_000, 10_000_000)
    for old, new in [
        ("min_ratio = 0.5", "min_ratio = 1.6"),
        ("background_factor = 3.0", "background_factor = 0.5"),
    ]:
        with pytest.raises(ValueError, match=r"\[hhv6\]"):
            load_grading(_variant(tmp_path, old, new))
    with pytest.raises(ValueError, match=r"\[roh\]"):
        load_grading(
            _variant(tmp_path, "isodisomy_min = 10000000", "isodisomy_min = 1")
        )
    with pytest.raises(ValueError, match=r"missing sections.*'roh'"):
        load_grading(_variant(tmp_path, "[roh]", "[unused_roh]"))


def _fake_checkout(path: Path) -> Path:
    (path / "src" / "genokeeper").mkdir(parents=True)
    (path / "pyproject.toml").write_text("")
    return path


def test_checkout_keeps_repo_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkout = _fake_checkout(tmp_path / "checkout")
    monkeypatch.setattr(config, "_REPO", checkout)
    cwd = tmp_path / "elsewhere"
    cwd.mkdir()
    assert load_paths(cwd=cwd).repo == checkout


def test_outside_checkout_uses_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "site-packages"
    empty.mkdir()
    monkeypatch.setattr(config, "_REPO", empty)
    paths = load_paths(cwd=tmp_path)
    assert paths.reports == tmp_path / "reports"


def test_demo_folder_lab_wins_over_home(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / DEMO_MARKER).write_text("{}")
    monkeypatch.setenv("GENOKEEPER_HOME", "/elsewhere")
    assert load_paths(tmp_path).lab == tmp_path / "lab"


def test_home_moves_lab_outside_demo_folder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GENOKEEPER_HOME", "/elsewhere")
    assert load_paths(tmp_path).lab == Path("/elsewhere")
