"""A miniature Genomics England PanelApp panel, shaped as the API writes panel 285: nothing
but what `genokeeper.sources.panelapp.load_panel` reads, under invented gene rows."""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from genokeeper.sources.panelapp import PANEL_ID, PANEL_VERSION, panelapp_path

# gene symbol, confidence_level, mode_of_inheritance, Ensembl gene id
type GeneRow = tuple[str, str, str, str]
AR = "BIALLELIC, autosomal or pseudoautosomal"
GENES: tuple[GeneRow, ...] = (
    ("BRCA1", "3", AR, "ENSG00000012048"),
    ("GENEGREEN", "3", AR, "ENSG01"),
    ("GENEAMBER", "2", "Unknown", "ENSG02"),
    ("GENERED", "1", "", "ENSG03"),
    ("GENENONE", "0", AR, "ENSG04"),
)


def panel_document(genes: Sequence[GeneRow] = GENES, **fields: Any) -> dict[str, Any]:
    """The panel's JSON document; `fields` replace its top-level keys."""
    return {
        "id": PANEL_ID,
        "name": "Intellectual disability",
        "version": PANEL_VERSION,
        "genes": [
            {
                "gene_data": {
                    "gene_symbol": symbol,
                    "ensembl_genes": {"GRch38": {"90": {"ensembl_id": gene_id}}},
                },
                "confidence_level": confidence,
                "mode_of_inheritance": inheritance,
            }
            for symbol, confidence, inheritance, gene_id in genes
        ],
    } | fields


def write_panel(path: Path, genes: Sequence[GeneRow] = GENES, **fields: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(panel_document(genes, **fields)))
    return path


def write_lab_panel(lab: Path, genes: Sequence[GeneRow] = GENES) -> Path:
    """The panel file where the bootstrap puts it under `lab`."""
    return write_panel(panelapp_path(lab), genes)
