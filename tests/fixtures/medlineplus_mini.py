"""MedlinePlus Genetics summaries in the published shape of `ghr-summaries.xml`: a default
namespace, the `html:` paragraphs, and names nested in related-gene and related-condition
lists that a loader must not take for the summary's own."""

from pathlib import Path

_OPEN = (
    '<?xml version="1.0" encoding="utf-8"?>\n'
    '<summaries xmlns="https://medlineplus.gov/download/ghr-summaries-20250602.xsd"'
    ' xmlns:html="http://www.w3.org/1999/xhtml">\n'
)
_PAGE = "https://medlineplus.gov/genetics"


def _html(paragraphs: tuple[str, ...]) -> str:
    return "".join(f"<html:p>{p}</html:p>" for p in paragraphs)


def gene(symbol: str, *paragraphs: str, published: str = "2025-08-25") -> str:
    """A gene summary whose `function` text is `paragraphs` (inline markup allowed)."""
    return (
        f"<gene-summary id='1'><gene-symbol>{symbol}</gene-symbol>"
        f"<name>{symbol.lower()} protein</name>"
        f"<ghr-page>{_PAGE}/gene/{symbol.lower()}</ghr-page>"
        "<text-list><text><text-role>function</text-role>"
        f"<html>{_html(paragraphs)}</html></text></text-list>"
        "<related-health-condition-list><related-health-condition>"
        f"<name>{symbol} related disorder</name>"
        "</related-health-condition></related-health-condition-list>"
        f"<synonym-list><synonym>{symbol} disorder</synonym></synonym-list>"
        f"<reviewed>2025-08</reviewed><published>{published}</published>"
        "</gene-summary>\n"
    )


def condition(
    page: str,
    name: str,
    *paragraphs: str,
    synonyms: tuple[str, ...] = (),
    published: str = "2020-09-08",
) -> str:
    """A condition summary at `<ghr-page>/condition/<page>` whose `description` text is
    `paragraphs`."""
    listed = "".join(f"<synonym>{s}</synonym>" for s in synonyms)
    return (
        f"<health-condition-summary id='2'><name>{name}</name>"
        f"<ghr-page>{_PAGE}/condition/{page}</ghr-page>"
        "<text-list><text><text-role>description</text-role>"
        f"<html>{_html(paragraphs)}</html></text></text-list>"
        "<related-gene-list><related-gene><gene-symbol>RELATED1</gene-symbol>"
        "<name>a related gene</name></related-gene></related-gene-list>"
        f"<synonym-list>{listed}</synonym-list>"
        f"<reviewed>2019-09</reviewed><published>{published}</published>"
        "</health-condition-summary>\n"
    )


def write_medlineplus(path: Path, *summaries: str) -> Path:
    """`summaries` as a `ghr-summaries.xml` at `path`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_OPEN + "".join(summaries) + "</summaries>\n", encoding="utf-8")
    return path
