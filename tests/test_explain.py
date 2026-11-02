"""Generated explanations: grounded prompts, keyed by model, prompt version and the prompt,
answers kept only when they cite the prompt's sources. The model is never called: every test
answers at the `Ask` seam, and `claude_ask`'s tests replace `subprocess.run`."""

import json
import os
import re
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures.medlineplus_mini import condition, gene, write_medlineplus
from fixtures.store_mini import Row, seed

from genokeeper import explain as explainer
from genokeeper.cli import main
from genokeeper.evidence.medlineplus import (
    MEDLINEPLUS_FILE,
    MedlinePlus,
    load_medlineplus,
)
from genokeeper.explain import (
    JOBS,
    MODEL,
    ExplainError,
    Explanation,
    Rejected,
    Source,
    cache_key,
    claude_ask,
    claude_command,
    explain,
    grounding,
    lookup,
    parse,
    prompt,
)
from genokeeper.lang import LANGUAGES
from genokeeper.query import HAPLOID_DOSAGE_NOTE, Filters, query
from genokeeper.store.db import Level, Validity, connect

GENE_CALL = {
    "diplotype": "*1/*2",
    "phenotype": "Intermediate Metabolizer",
    "activity_score": 1.5,
    "callers": {"BCyrius": "*1/*2", "Aldy": "*1/*2"},
    "agreed": True,
    "level": "high",
    "warnings": ["a gene call warning"],
}
ROWS = [
    Row(
        "clinvar:1",
        gene="HFE",
        condition="Hereditary_hemochromatosis|not_provided",
        warnings=["a call warning"],
    ),
    Row("clinvar:2", gene="BRCA1", pos=1002),
    Row(
        "pharmcat:CYP2D6",
        gene="CYP2D6",
        level=None,
        topic="drugs",
        gene_call=GENE_CALL,
        claim={"schema": "genomelab.claim/1", "type": "pharmcat_gene", "drugs": []},
    ),
]
HFE_TEXT = "The HFE gene makes a protein."
HH_TEXT = "Hereditary hemochromatosis is a disorder."


def _medlineplus(path: Path, published: str = "2020-09-08") -> Path:
    return write_medlineplus(
        path,
        gene("HFE", HFE_TEXT),
        condition(
            "hereditary-hemochromatosis",
            "Hereditary hemochromatosis",
            HH_TEXT,
            published=published,
        ),
    )


@pytest.fixture
def store(tmp_path: Path) -> Path:
    path = tmp_path / "genomelab.duckdb"
    with connect(path) as con:
        seed(con, ROWS)
    return path


def _findings(store: Path) -> dict[str, dict[str, Any]]:
    with duckdb.connect(str(store), read_only=True) as con:
        loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
        return {f["id"]: f for f in query(con, None, loosest).findings}


@pytest.fixture
def findings(store: Path) -> dict[str, dict[str, Any]]:
    return _findings(store)


@pytest.fixture
def medlineplus(tmp_path: Path) -> MedlinePlus:
    return load_medlineplus(_medlineplus(tmp_path / "ghr.xml"))


def _document(text: str) -> dict[str, Any]:
    """A prompt's facts and sources, the JSON after the instructions."""
    return json.loads(text[text.index("\n{\n") + 1 :])


class Recorder:
    """An `Ask` that records each prompt and answers with `answer(prompt)`."""

    def __init__(self, answer: Any = None) -> None:
        self.prompts: list[str] = []
        self.answer = answer or (
            lambda p: json.dumps(
                {
                    "title": "Заголовок.",
                    "text": "Объяснение.",
                    "cites": [_document(p)["sources"][0]["id"]],
                },
                ensure_ascii=False,
            )
        )

    def __call__(self, text: str) -> str:
        self.prompts.append(text)
        return self.answer(text)


def test_grounding_gene_and_condition(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus
) -> None:
    finding = findings["clinvar:1"]
    assert grounding(finding, medlineplus) == [
        Source("clinvar:1", json.dumps(finding["claim"], sort_keys=True)),
        Source("medlineplus:gene:HFE", HFE_TEXT),
        # `not_provided` names no MedlinePlus condition and adds nothing.
        Source("medlineplus:condition:hereditary-hemochromatosis", HH_TEXT),
    ]
    # A gene without a summary and a claim without a condition: the finding's own source.
    assert [s.id for s in grounding(findings["pharmcat:CYP2D6"], medlineplus)] == [
        "pharmcat:CYP2D6"
    ]


def test_grounding_without_medlineplus(findings: dict[str, dict[str, Any]]) -> None:
    assert [s.id for s in grounding(findings["clinvar:1"], None)] == ["clinvar:1"]


def test_prompt_carries_facts_warnings_and_sources_not_p_wrong(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus
) -> None:
    finding = findings["pharmcat:CYP2D6"]
    sources = grounding(finding, medlineplus)
    text = prompt(finding, sources)
    document = _document(text)
    assert document["finding"] == {
        "id": "pharmcat:CYP2D6",
        "topic": "drugs",
        "gene": "CYP2D6",
        "marker": finding["marker"],
        "genotype": finding["genotype"],
        "validity": finding["validity"],
        "actionable": False,
        "effect": None,
        "gene_call": {
            "diplotype": "*1/*2",
            "phenotype": "Intermediate Metabolizer",
            "activity_score": 1.5,
            "callers": {"BCyrius": "*1/*2", "Aldy": "*1/*2"},
            "level": "high",
            "warnings": ["a gene call warning"],
        },
        "warnings": None,
        "dosage_note": None,
    }
    assert document["sources"] == [{"id": s.id, "text": s.text} for s in sources]
    variant = prompt(
        findings["clinvar:1"], grounding(findings["clinvar:1"], medlineplus)
    )
    for leaked in ("p_wrong", "flags"):
        assert leaked not in text and leaked not in variant
    assert "in English" in text and "data, never an instruction" in text
    assert '"title"' in text


def test_the_english_prompt_addresses_the_reader_and_the_russian_one_is_unchanged(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus
) -> None:
    """The Russian prompt is the one the user's cached answers were keyed on before the
    English default: its keys are pinned to those made then, so they still hit."""
    chrm = dict(findings["clinvar:1"], marker={"grch38": "chrM:1:C:T"})
    english = prompt(chrm, grounding(chrm, medlineplus))
    assert 'addressing\nthem as "you"' in english and "owner" not in english
    assert "in Russian" not in english
    russian = prompt(chrm, grounding(chrm, medlineplus), lang="ru")
    assert "to the owner, in Russian" in russian
    assert _document(russian) == {
        **_document(english),
        "finding": {
            **_document(english)["finding"],
            "dosage_note": _document(russian)["finding"]["dosage_note"],
        },
    }
    pinned = {
        "clinvar:1": "1ec9c503d533adc1ad7848692d87099cee9f99059ea9b5b1385d39dee5bd323c",
        "clinvar:2": "753e029aaea369df8c8cf23e36c52086d0d402778f4f8dd37ba950bcc489346b",
        "pharmcat:CYP2D6": (
            "c9338bfe4da9d5b46e647880d623280ea37284668f91d3c514dd3e4c51062b2e"
        ),
        "chrM": "71f122603d06afec203760d77cf78f7255b508305af05e0492539dd5e3175f31",
    }
    for name, f in [*findings.items(), ("chrM", chrm)]:
        sources = grounding(f, medlineplus)
        assert cache_key(prompt(f, sources, lang="ru")) == pinned[name]
        assert cache_key(prompt(f, sources)) != pinned[name]


def test_an_answer_is_never_served_in_the_other_language(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus, tmp_path: Path
) -> None:
    cache = tmp_path / "explanations"
    russian = Recorder()
    counts = explain(list(findings.values()), medlineplus, cache, russian, lang="ru")
    assert counts["written"] == 3
    finding = findings["clinvar:1"]
    assert lookup(finding, medlineplus, cache, lang="ru") is not None
    assert lookup(finding, medlineplus, cache) is None
    english = Recorder(
        lambda p: json.dumps(
            {
                "title": "Title.",
                "text": "Explanation.",
                "cites": [_document(p)["sources"][0]["id"]],
            }
        )
    )
    counts = explain(list(findings.values()), medlineplus, cache, english)
    assert counts == {"cached": 0, "written": 3, "failed": 0}
    assert all("in English" in text for text in english.prompts)
    found = lookup(finding, medlineplus, cache)
    assert found is not None and found.title == "Title."
    again = lookup(finding, medlineplus, cache, lang="ru")
    assert again is not None and again.title == "Заголовок."


def test_a_gene_call_warning_changes_the_key(
    findings: dict[str, dict[str, Any]],
) -> None:
    """A caveat on the gene call (a repeat allele beyond the fragment length) reaches the
    model, so its appearing makes a new key."""
    f = findings["pharmcat:CYP2D6"]
    caveat = "allele beyond the fragment length"
    warned = dict(f, gene_call={**f["gene_call"], "warnings": [caveat]})
    assert caveat in prompt(warned, grounding(warned, None))
    assert cache_key(prompt(warned, grounding(warned, None))) != cache_key(
        prompt(f, grounding(f, None))
    )


def test_no_absolute_local_path_reaches_the_prompt(
    findings: dict[str, dict[str, Any]],
) -> None:
    """A repeat claim's REViewer image path names the user's home and checkout: the
    prompt (sent to Anthropic) and its key leave it out."""
    f = findings["pharmcat:CYP2D6"]
    image = "/Users/someone/dna/reports/work/repeats/k/rev.HTT.svg"
    claim = {"type": "repeat_expansion", "gene": "HTT", "motif": "CAG"}
    imaged = dict(f, claim={**claim, "image": image})
    text = prompt(imaged, grounding(imaged, None))
    assert "/Users/" not in text and "rev.HTT.svg" not in text and "CAG" in text
    plain = dict(f, claim=claim)
    assert cache_key(text) == cache_key(prompt(plain, grounding(plain, None)))


def test_prompt_ignores_stored_key_order(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus
) -> None:
    finding = findings["clinvar:1"]
    sources = grounding(finding, medlineplus)
    reordered = {
        k: dict(reversed(v.items())) if isinstance(v, dict) else v
        for k, v in finding.items()
    }
    assert prompt(reordered, sources) == prompt(finding, sources)


def test_prompt_says_what_dosage_means_off_the_autosomes(
    findings: dict[str, dict[str, Any]],
) -> None:
    def facts(contig: str, warnings: list[str]) -> dict[str, Any]:
        f = dict(findings["clinvar:1"], marker={"grch38": f"{contig}:1:C:T"})
        f["reliability"] = {**f["reliability"], "warnings": warnings}
        return _document(prompt(f, grounding(f, None)))["finding"]

    assert facts("chr17", [])["dosage_note"] is None
    assert facts("chrX", [])["dosage_note"] == HAPLOID_DOSAGE_NOTE
    chrm = facts("chrM", ["heteroplasmy at level 0.23"])
    assert "present in the sample's mitochondria" in chrm["dosage_note"]
    assert chrm["warnings"] == ["heteroplasmy at level 0.23"]
    f = dict(findings["clinvar:1"], marker={"grch38": "chrM:1:C:T"})
    russian = _document(prompt(f, grounding(f, None), lang="ru"))["finding"]
    assert "present in the owner's mitochondria" in russian["dosage_note"]


def _with_version(finding: dict[str, Any], version: str) -> dict[str, Any]:
    """`finding` as a newer release states it: its own source's version changed."""
    return dict(finding, sources=[dict(finding["sources"][0], version=version)])


def test_a_new_release_alone_keeps_the_key(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus, tmp_path: Path
) -> None:
    """A release changes an explanation only when the finding's facts or its sources' text
    changed: a new version of the finding's own source, or a republished MedlinePlus
    summary with the same text, is the same key."""
    f = findings["clinvar:1"]
    key = cache_key(prompt(f, grounding(f, medlineplus)))
    newer = _with_version(f, "20261123")
    assert cache_key(prompt(newer, grounding(newer, medlineplus))) == key
    republished = load_medlineplus(_medlineplus(tmp_path / "new.xml", "2026-01-01"))
    assert cache_key(prompt(f, grounding(f, republished))) == key


def test_changed_text_changes_the_key(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus, tmp_path: Path
) -> None:
    """The claim's text and a summary's text are what the answer is grounded in."""
    f = findings["clinvar:1"]
    key = cache_key(prompt(f, grounding(f, medlineplus)))
    claim = dict(f["claim"], condition="Hereditary_hemochromatosis|something_else")
    reclaimed = dict(f, claim=claim)
    assert cache_key(prompt(reclaimed, grounding(reclaimed, medlineplus))) != key
    summary = load_medlineplus(
        write_medlineplus(
            tmp_path / "text.xml",
            gene("HFE", HFE_TEXT + " It was revised."),
            condition(
                "hereditary-hemochromatosis",
                "Hereditary hemochromatosis",
                HH_TEXT,
            ),
        )
    )
    assert cache_key(prompt(f, grounding(f, summary))) != key


def test_key_changes_with_model_and_prompt_version(
    findings: dict[str, dict[str, Any]],
    medlineplus: MedlinePlus,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    finding = findings["clinvar:1"]
    text = prompt(finding, grounding(finding, medlineplus))
    key = cache_key(text)
    assert key == cache_key(prompt(finding, grounding(finding, medlineplus)))
    monkeypatch.setattr(explainer, "MODEL", "claude-other-1")
    assert cache_key(text) != key
    monkeypatch.undo()
    monkeypatch.setattr(explainer, "PROMPT_VERSION", explainer.PROMPT_VERSION + 1)
    assert cache_key(text) != key


def test_the_prompt_version_is_three() -> None:
    """The prompt no longer carries a source's version, so every key of version 2 is
    stale by construction."""
    assert explainer.PROMPT_VERSION == 3


ALLOWED = {"clinvar:1", "medlineplus:gene:HFE"}


def test_parse_keeps_a_cited_answer() -> None:
    answer = {"title": "Заголовок.", "text": "Текст.", "cites": ["clinvar:1"]}
    assert parse(json.dumps(answer), ALLOWED) == Explanation(
        "Заголовок.", "Текст.", ["clinvar:1"], MODEL
    )


@pytest.mark.parametrize("title", [None, "", " ", 3])
def test_rejects_no_title(title: Any) -> None:
    answer: dict[str, Any] = {"text": "Текст.", "cites": ["clinvar:1"]}
    if title is not None:
        answer["title"] = title
    assert parse(json.dumps(answer), ALLOWED) == Rejected("no_title")


@pytest.mark.parametrize(
    "answer",
    [
        '{"text": "Текст.", "cites": []}',
        '{"text": "Текст."}',
        '{"text": "Текст.", "cites": "clinvar:1"}',
    ],
)
def test_rejects_no_cites(answer: str) -> None:
    assert parse(answer, ALLOWED) == Rejected("no_cites")


@pytest.mark.parametrize(
    "cites", [["clinvar:2"], ["clinvar:1", "medlineplus:gene:BRCA1"], [{"id": 1}]]
)
def test_rejects_unknown_cite(cites: list[Any]) -> None:
    answer = json.dumps({"text": "Текст.", "cites": cites})
    assert parse(answer, ALLOWED) == Rejected("unknown_cite")


@pytest.mark.parametrize(
    ("answer", "reason"),
    [
        # Only one ```json fence around the whole answer is dropped.
        (
            '```json\n{"text": "Текст.", "cites": ["clinvar:1"]}\n```\nПояснение.',
            "not_json",
        ),
        ('```\n{"text": "Текст.", "cites": ["clinvar:1"]}\n```', "not_json"),
        ('{"text": "Текст.", "cites": ["clinvar:1"]} and more', "not_json"),
        ('["clinvar:1"]', "not_json"),
        ('{"text": "", "cites": ["clinvar:1"]}', "no_text"),
        ('{"text": 3, "cites": ["clinvar:1"]}', "no_text"),
    ],
)
def test_rejects_not_json(answer: str, reason: str) -> None:
    assert parse(answer, ALLOWED) == Rejected(reason)


def test_rejected_answers_counted_nothing_written(
    findings: dict[str, dict[str, Any]],
    medlineplus: MedlinePlus,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cache = tmp_path / "explanations"
    sent = [
        '{"text": "Текст.", "cites": []}',
        '{"text": "Текст.", "cites": ["clinvar:999"]}',
        "Не JSON.\nВторая строка.",
    ]
    answers = iter(sent)
    counts = explain(
        list(findings.values()),
        medlineplus,
        cache,
        Recorder(lambda _: next(answers)),
        jobs=1,
    )
    # One stderr line per finding asked, a rejected answer's text repr'd onto it.
    reasons = ("no_cites", "unknown_cite", "not_json")
    assert capsys.readouterr().err.splitlines() == [
        f"genokeeper explain: {i}/3 {id_} rejected:{reason} {answer!r}"
        for i, (id_, reason, answer) in enumerate(zip(findings, reasons, sent), 1)
    ]
    assert counts == {
        "cached": 0,
        "written": 0,
        "failed": 0,
        "rejected:no_cites": 1,
        "rejected:unknown_cite": 1,
        "rejected:not_json": 1,
    }
    assert not cache.exists() or not list(cache.iterdir())
    untitled = Recorder(
        lambda p: json.dumps(
            {"text": "Т.", "cites": [_document(p)["sources"][0]["id"]]}
        )
    )
    counts = explain(list(findings.values()), medlineplus, cache, untitled)
    assert counts["rejected:no_title"] == 3 and counts["written"] == 0
    assert not cache.exists() or not list(cache.iterdir())


def test_cached_answer_not_asked_again(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus, tmp_path: Path
) -> None:
    cache = tmp_path / "explanations"
    first = Recorder()
    counts = explain(list(findings.values()), medlineplus, cache, first)
    assert counts == {"cached": 0, "written": 3, "failed": 0}
    finding = findings["clinvar:1"]
    assert lookup(finding, medlineplus, cache) == Explanation(
        "Заголовок.", "Объяснение.", ["clinvar:1"], MODEL
    )
    again = Recorder()
    assert explain(list(findings.values()), medlineplus, cache, again) == {
        "cached": 3,
        "written": 0,
        "failed": 0,
    }
    assert again.prompts == []


def test_failed_call_counted_and_run_goes_on(
    findings: dict[str, dict[str, Any]],
    medlineplus: MedlinePlus,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cache = tmp_path / "explanations"
    calls = Recorder()

    def flaky(text: str) -> str:
        if not calls.prompts:
            calls.prompts.append(text)
            raise ExplainError("claude exited 1")
        return calls(text)

    counts = explain(list(findings.values()), medlineplus, cache, flaky, jobs=1)
    assert counts == {"cached": 0, "written": 2, "failed": 1}
    assert len(list(cache.iterdir())) == 2
    first, *rest = findings
    assert capsys.readouterr().err.splitlines() == [
        f"genokeeper explain: 1/3 {first} failed: 'claude exited 1'",
        *(f"genokeeper explain: {i}/3 {id_} written" for i, id_ in enumerate(rest, 2)),
    ]
    # A rerun asks only the failed one again; cached findings print nothing.
    explain(list(findings.values()), medlineplus, cache, Recorder())
    assert capsys.readouterr().err == f"genokeeper explain: 1/1 {first} written\n"


def test_a_multi_line_failure_reason_stays_on_its_line(
    findings: dict[str, dict[str, Any]],
    medlineplus: MedlinePlus,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def failing(_: str) -> str:
        raise ExplainError("claude exited 1: Traceback\n  line 2")

    explain(list(findings.values())[:1], medlineplus, tmp_path, failing)
    reason = repr("claude exited 1: Traceback\n  line 2")
    assert capsys.readouterr().err.splitlines() == [
        f"genokeeper explain: 1/1 {next(iter(findings))} failed: {reason}"
    ]


def _own_id(text: str) -> str:
    return _document(text)["finding"]["id"]


def _many(tmp_path: Path, n: int) -> list[dict[str, Any]]:
    path = tmp_path / "many.duckdb"
    with connect(path) as con:
        seed(con, [Row(f"clinvar:{i}", pos=1000 + i) for i in range(1, n + 1)])
    return list(_findings(path).values())


def test_jobs_bounds_the_calls_at_once(
    medlineplus: MedlinePlus, tmp_path: Path
) -> None:
    # Each call waits until `jobs` calls have run at once, so a pool narrower than `jobs`
    # never gets there (caught: `max_workers=1`) and a wider one exceeds it (caught:
    # `max_workers` not `jobs`).
    jobs, running, peak = 2, 0, 0
    changed = threading.Condition()

    def ask(text: str) -> str:
        nonlocal running, peak
        with changed:
            running += 1
            peak = max(peak, running)
            changed.notify_all()
            changed.wait_for(lambda: peak >= jobs, timeout=5)
        time.sleep(0.05)
        with changed:
            running -= 1
        return Recorder().answer(text)

    counts = explain(_many(tmp_path, 6), medlineplus, tmp_path / "e", ask, jobs=jobs)
    assert counts["written"] == 6
    assert peak == jobs


def test_parallel_answers_all_cached_and_counted(
    medlineplus: MedlinePlus, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Outcomes keyed by the finding, so the order the pool runs them in cannot matter.
    # Caught: an answer not written or written under another finding's key, an outcome
    # counted under another name, a line lost or numbered twice.
    findings = _many(tmp_path, 12)
    outcome = {
        f["id"]: ("failed", "rejected:not_json", "written", "written")[i % 4]
        for i, f in enumerate(findings)
    }

    def ask(text: str) -> str:
        time.sleep(0.01)
        if outcome[_own_id(text)] == "failed":
            raise ExplainError("claude exited 1: API error 429 usage_limit_reached")
        if outcome[_own_id(text)] == "rejected:not_json":
            return "Не JSON."
        return Recorder().answer(text)

    cache = tmp_path / "explanations"
    counts = explain(findings, medlineplus, cache, ask, jobs=JOBS)
    assert counts == {
        "cached": 0,
        "written": 6,
        "failed": 3,
        "rejected:not_json": 3,
    }
    for f in findings:
        cached = lookup(f, medlineplus, cache)
        assert (cached is not None) == (outcome[f["id"]] == "written")
        assert cached is None or cached.cites == [f["id"]]
    assert sorted(p.suffix for p in cache.iterdir()) == [".json"] * 6
    lines = capsys.readouterr().err.splitlines()
    seen = [
        re.fullmatch(r"genokeeper explain: (\d+)/12 (\S+) (\S+).*", x) for x in lines
    ]
    assert all(seen)
    assert sorted(int(m[1]) for m in seen if m) == list(range(1, 13))
    assert {m[2]: m[3] for m in seen if m} == {
        id_: f"{o}:" if o == "failed" else o for id_, o in outcome.items()
    }
    # A rerun asks only the ones not written.
    again = Recorder()
    counts = explain(findings, medlineplus, cache, again, jobs=JOBS)
    assert counts == {"cached": 6, "written": 6, "failed": 0}
    assert sorted(map(_own_id, again.prompts)) == sorted(
        id_ for id_, o in outcome.items() if o != "written"
    )


def test_an_error_in_one_call_keeps_the_running_calls_answers(
    findings: dict[str, dict[str, Any]], medlineplus: MedlinePlus, tmp_path: Path
) -> None:
    # The first call raises an error that is not `ExplainError` while the second is still
    # waiting on its answer; that answer is cached before the error reaches the caller.
    # Caught: a pool shut down without waiting for its running calls.
    first, second, _ = findings.values()
    started, raised = threading.Event(), threading.Event()

    def ask(text: str) -> str:
        if _own_id(text) == first["id"]:
            assert started.wait(timeout=5)
            raised.set()
            raise RuntimeError("a bug in one call")
        started.set()
        assert raised.wait(timeout=5)
        time.sleep(0.2)
        return Recorder().answer(text)

    cache = tmp_path / "explanations"
    with pytest.raises(RuntimeError, match="a bug in one call"):
        explain([first, second], medlineplus, cache, ask, jobs=2)
    assert lookup(second, medlineplus, cache) is not None
    assert lookup(first, medlineplus, cache) is None


def test_claude_command() -> None:
    assert claude_command() == [
        "claude",
        "-p",
        "--model",
        "claude-opus-5-5",
        "--tools",
        "",
        "--strict-mcp-config",
        "--mcp-config",
        '{"mcpServers":{}}',
        "--setting-sources",
        "",
        "--disable-slash-commands",
        "--no-session-persistence",
        "--output-format",
        "json",
    ]
    assert "--bare" not in claude_command()


# Claude Code's `--output-format json` result, as its documentation gives it.
RESULT = {
    "type": "result",
    "subtype": "success",
    "is_error": False,
    "duration_ms": 2300,
    "num_turns": 1,
    "result": '{"text": "Текст.", "cites": ["clinvar:1"]}',
    "session_id": "00000000-0000-0000-0000-000000000000",
    "total_cost_usd": 0.01,
}


def _fake_run(
    monkeypatch: pytest.MonkeyPatch, returncode: int, stdout: str
) -> list[dict[str, Any]]:
    calls: list[dict[str, Any]] = []

    def run(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append({"argv": argv, **kwargs, "listed": os.listdir(kwargs["cwd"])})
        return subprocess.CompletedProcess(argv, returncode, stdout, "an error line\n")

    monkeypatch.setattr(subprocess, "run", run)
    return calls


def test_claude_ask_reads_result(monkeypatch: pytest.MonkeyPatch) -> None:
    leaking = ("CLAUDECODE", "ANTHROPIC_BASE_URL", "ANTHROPIC_API_KEY")
    for name in leaking:
        monkeypatch.setenv(name, "from the parent")
    calls = _fake_run(monkeypatch, 0, json.dumps(RESULT))
    assert claude_ask("the prompt") == RESULT["result"]
    [call] = calls
    assert call["argv"] == claude_command()
    assert call["input"] == "the prompt"  # on stdin, never in the argv
    assert call["listed"] == [] and not Path(call["cwd"]).exists()
    assert call["timeout"] > 0
    assert "HOME" in call["env"] and not set(leaking) & set(call["env"])


@pytest.mark.parametrize(
    ("returncode", "stdout"),
    [
        (
            0,
            json.dumps(
                {**RESULT, "is_error": True, "subtype": "error_during_execution"}
            ),
        ),
        (1, json.dumps(RESULT)),
        (0, "not json"),
        (0, json.dumps({**RESULT, "result": None})),
    ],
)
def test_claude_ask_refuses_an_error(
    monkeypatch: pytest.MonkeyPatch, returncode: int, stdout: str
) -> None:
    _fake_run(monkeypatch, returncode, stdout)
    with pytest.raises(ExplainError):
        claude_ask("the prompt")


def test_claude_ask_names_the_result_of_a_failed_process(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Claude Code 2.1 exits 1 on an expired login with an empty stderr and the reason in
    # its JSON result (a run printed "claude exited 1: " for every finding).
    reason = "Failed to authenticate: OAuth session expired and could not be refreshed"
    _fake_run(
        monkeypatch, 1, json.dumps({**RESULT, "is_error": True, "result": reason})
    )
    with pytest.raises(ExplainError, match=f"claude exited 1: {reason}"):
        claude_ask("the prompt")


@pytest.mark.parametrize(
    "result",
    [RESULT, {**RESULT, "is_error": True, "result": ""}],
    ids=["an answer, not an error", "an empty reason"],
)
def test_a_failed_process_without_an_error_reason_names_its_stderr(
    monkeypatch: pytest.MonkeyPatch, result: dict[str, Any]
) -> None:
    _fake_run(monkeypatch, 1, json.dumps(result))
    with pytest.raises(ExplainError, match="^claude exited 1: an error line$"):
        claude_ask("the prompt")


@pytest.mark.parametrize(
    ("returncode", "result", "reason"),
    [
        # A usage limit on 2.1.304 (Claude Code issue 103807, 2026-10-20): no `result`.
        (
            1,
            {
                "type": "result",
                "subtype": "success",
                "terminal_reason": "api_error",
                "is_error": True,
                "api_error_status": 429,
                "api_error": "usage_limit_reached",
            },
            "claude exited 1: API error 429 usage_limit_reached",
        ),
        # A rate limit that exits 0 (issue 79500, still reported on 2.1.267, 2026-09-15).
        (
            0,
            {**RESULT, "is_error": True, "result": "API Error: Rate limit reached"},
            "claude answered with an error: API Error: Rate limit reached",
        ),
        (
            0,
            {
                "type": "result",
                "subtype": "success",
                "is_error": True,
                "api_error_status": 529,
                "api_error": "overloaded",
            },
            "claude answered with an error: API error 529 overloaded",
        ),
    ],
    ids=["usage limit, exit 1", "rate limit, exit 0", "overloaded, no text"],
)
def test_a_rate_limit_is_a_failed_call_with_its_reason(
    monkeypatch: pytest.MonkeyPatch, returncode: int, result: Any, reason: str
) -> None:
    # Caught: the API error fields not read (the reason falls back to stderr or stdout).
    _fake_run(monkeypatch, returncode, json.dumps(result))
    with pytest.raises(ExplainError) as raised:
        claude_ask("the prompt")
    assert str(raised.value) == reason


def _lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    lab = tmp_path / "lab"
    lab.mkdir()
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    return lab


def test_dry_run_calls_nothing_and_prints_prompts(
    store: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    write_medlineplus(
        lab / "medlineplus" / "20261001" / MEDLINEPLUS_FILE,
        gene("HFE", "An older summary."),
    )
    _medlineplus(lab / "medlineplus" / "20261020" / MEDLINEPLUS_FILE)

    def refuse(*_: Any, **__: Any) -> Any:
        raise AssertionError("a dry run called out")

    monkeypatch.setattr(explainer, "claude_ask", refuse)
    monkeypatch.setattr(subprocess, "run", refuse)
    cache = tmp_path / "explanations"
    command = ["explain", "--db", str(store), "--cache", str(cache)]
    assert main([*command, "--dry-run"]) == 0
    printed = capsys.readouterr().out
    assert not cache.exists()
    # A real run, answered by a stub, sends exactly the prompts the dry run printed.
    sent = Recorder()
    monkeypatch.setattr(explainer, "claude_ask", sent)
    assert main(command) == 0
    assert len(sent.prompts) == 3
    for text in sent.prompts:
        assert text in printed
    # The newest dated folder's summaries ground the prompts.
    assert HFE_TEXT in printed and "An older summary." not in printed
    assert printed.endswith(f"3 calls to {MODEL}\n")
    run = capsys.readouterr()
    assert run.out.splitlines() == ["cached: 0", "written: 3", "failed: 0"]
    assert len(run.err.splitlines()) == 3
    monkeypatch.setattr(explainer, "claude_ask", refuse)
    assert main([*command, "--dry-run"]) == 0
    assert capsys.readouterr().out == f"0 calls to {MODEL}\n"


def test_dry_run_prints_the_prompts_of_the_chosen_language(
    store: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    _medlineplus(lab / "medlineplus" / "20261020" / MEDLINEPLUS_FILE)
    command = ["explain", "--db", str(store), "--cache", str(tmp_path / "e")]

    def printed(*argv: str) -> str:
        assert main([*command, "--dry-run", *argv]) == 0
        return capsys.readouterr().out

    monkeypatch.delenv("GENOKEEPER_LANG", raising=False)
    assert "in English" in printed() and "in Russian" not in printed()
    assert "in Russian" in printed("--lang", "ru")
    monkeypatch.setenv("GENOKEEPER_LANG", "ru")
    assert "in Russian" in printed()
    assert "in English" in printed("--lang", "en")  # the option wins
    # A real run asks in the chosen language too.
    langs: list[str] = []

    def record(*_: Any, lang: str, **__: Any) -> dict[str, int]:
        langs.append(lang)
        return {"cached": 0, "written": 0, "failed": 0}

    monkeypatch.setattr(explainer, "explain", record)
    assert main(command) == 0
    assert main([*command, "--lang", "en"]) == 0
    assert langs == ["ru", "en"]
    with pytest.raises(SystemExit):
        main([*command, "--lang", "de"])
    assert "invalid choice" in capsys.readouterr().err


def test_explain_exits_1_when_a_call_failed(
    store: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    _medlineplus(lab / "medlineplus" / "20261020" / MEDLINEPLUS_FILE)

    def fail(_: str) -> str:
        raise ExplainError("claude exited 1: Failed to authenticate")

    monkeypatch.setattr(explainer, "claude_ask", fail)
    cache = tmp_path / "explanations"
    assert main(["explain", "--db", str(store), "--cache", str(cache)]) == 1
    assert capsys.readouterr().out.splitlines() == [
        "cached: 0",
        "written: 0",
        "failed: 3",
    ]


def test_jobs_option_reaches_the_pool(
    store: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Caught: `--jobs` parsed and not passed on, or a default other than `JOBS`.
    lab = _lab(tmp_path, monkeypatch)
    _medlineplus(lab / "medlineplus" / "20261020" / MEDLINEPLUS_FILE)
    seen: list[int] = []

    def record(*_: Any, jobs: int, **__: Any) -> dict[str, int]:
        seen.append(jobs)
        return {"cached": 0, "written": 0, "failed": 0}

    monkeypatch.setattr(explainer, "explain", record)
    command = ["explain", "--db", str(store), "--cache", str(tmp_path / "e")]
    assert main(command) == 0
    assert main([*command, "--jobs", "1"]) == 0
    assert seen == [JOBS, 1]
    with pytest.raises(SystemExit):
        main([*command, "--jobs", "0"])
    assert "not a positive number" in capsys.readouterr().err


def test_explain_refuses_without_medlineplus(
    store: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _lab(tmp_path, monkeypatch)
    cache = tmp_path / "explanations"
    assert main(["explain", "--db", str(store), "--cache", str(cache)]) == 2
    err = capsys.readouterr().err
    assert err.startswith("genokeeper explain: ") and "scripts/bootstrap.sh" in err
    assert len(err.splitlines()) == 1
    missing = tmp_path / "nowhere.duckdb"
    assert main(["explain", "--db", str(missing), "--dry-run"]) == 2
    assert str(missing) in capsys.readouterr().err
    assert not cache.exists()


def test_parse_drops_one_json_fence() -> None:
    answer = {"title": "Заголовок", "text": "Текст.", "cites": ["clinvar:1"]}
    fenced = f"```json\n{json.dumps(answer, ensure_ascii=False)}\n```\n"
    assert parse(fenced, ALLOWED) == Explanation(model=MODEL, **answer)


WATCH_REASON = "screening result, needs confirmation: borderline, indeterminate"
STAGE3_ROWS = [
    Row(
        "expansionhunter:HTT",
        kind="watch",
        gene="HTT",
        level=None,
        claim={"schema": "genomelab.claim/1", "type": "repeat_expansion"},
        reason=WATCH_REASON,
        gene_call={"diplotype": "CAG 17/25", "phenotype": "borderline"},
    ),
    Row(
        "expansionhunter:DMPK",
        kind="non_finding",
        gene="DMPK",
        level=None,
        claim={"schema": "genomelab.claim/1", "type": "repeat_expansion"},
    ),
]


def test_explain_asks_for_a_stage3_watch_item_as_a_watch_item(
    store: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A stage-3 watch item reaches neither the feed nor `query`'s findings, so `explain`
    asks for it apart, told it is a watch item, never a finding, in either language; a
    stage-3 non-finding is not asked for (its call and stored reason say it), nor is a
    watch of another source or one at a variant."""
    other = {"schema": "genomelab.claim/1", "type": "x"}
    with connect(store) as con:
        seed(con, STAGE3_ROWS)
        seed(
            con,
            [
                Row("clinvar:3", kind="watch", gene="TTN", pos=1003),
                Row("hhv6:ici", kind="watch", gene="U94", level=None, claim=other),
                Row("expansionhunter:ATN1", kind="watch", gene="ATN1", pos=1004,
                    claim=other),
            ],
        )  # fmt: skip
        for prefix in ("expansionhunter", "hhv6"):
            con.execute(
                "UPDATE statement SET source = ? WHERE statement_id LIKE ?",
                [prefix, f"{prefix}:%"],
            )
    lab = _lab(tmp_path, monkeypatch)
    _medlineplus(lab / "medlineplus" / "20261020" / MEDLINEPLUS_FILE)
    for lang in LANGUAGES:
        sent = Recorder()
        monkeypatch.setattr(explainer, "claude_ask", sent)
        cache = str(tmp_path / f"e-{lang}")
        command = ["explain", "--db", str(store), "--cache", cache, "--lang", lang]
        assert main(command) == 0
        capsys.readouterr()
        asked = {_document(p)["finding"]["id"]: p for p in sent.prompts}
        assert set(asked) == {
            "clinvar:1",
            "clinvar:2",
            "pharmcat:CYP2D6",
            "expansionhunter:HTT",
        }
        watch = asked.pop("expansionhunter:HTT")
        assert watch.index(explainer.WATCH_NOTE) < watch.index(
            "The finding and its sources:"
        )
        assert _document(watch)["finding"]["kind"] == "watch"
        for text in asked.values():
            assert explainer.WATCH_NOTE not in text
            assert "kind" not in _document(text)["finding"]
