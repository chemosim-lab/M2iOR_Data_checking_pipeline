# pipeline/scripts/process_assay.py
"""Free-text Assay columns respelled to a per-column registry of canonical
values, so a new study or delivery cannot silently reintroduce a casing/
spacing variant of an equipment or construct name already fixed elsewhere
(e.g. "OC-725C Amplifier" for the registry's "OC-725C").

Modeled on `process_molecules.registry_spelling` / `process_sources.
normalize_doi_column`: only a case or spacing difference from a known
registry entry is auto-fixed. A cell whose normalized value has no registry
match is new data, not a formatting problem, and is left untouched - there
is no fixed allowed list, so legitimately new equipment or constructs are
never blocked. The registry itself (`assay_registry.json`) is seeded from
the corpus of already-curated studies by `scripts/tools/build_assay_registry.py`.
"""

import json
import re
from pathlib import Path

import pandas as pd

from scripts.columns import (
    ASSAY,
    EXPRESSION_CELL_TYPE,
    EXPRESSION_DRIVER,
    EXPRESSION_SYSTEM,
    ODOR_DELIVERY,
    RECORDING_SYSTEM,
    TYPE,
)
from scripts.report import Issue, emit

REGISTRY_PATH = Path(__file__).resolve().parent.parent / "assay_registry.json"

# scripts/columns.py L45-52: the Assay group's free-text columns (as opposed
# to Experimental technique, a closed vocabulary already checked by
# process_responses._validate_allowed_values).
ASSAY_FORMAT_COLUMNS: list[str] = [
    RECORDING_SYSTEM,
    TYPE,
    EXPRESSION_SYSTEM,
    EXPRESSION_CELL_TYPE,
    EXPRESSION_DRIVER,
    ODOR_DELIVERY,
]


def registry_key(value: str) -> str:
    """Normalized form used to match a cell against the registry: case and
    surrounding/repeated whitespace only, so two spellings differing in
    substance (not just formatting) never collide."""
    return re.sub(r"\s+", " ", value.strip().lower())


def load_registry() -> dict[str, dict[str, str]]:
    """Column name -> registry key -> canonical spelling."""
    if not REGISTRY_PATH.exists():
        return {}
    with REGISTRY_PATH.open(encoding="utf-8") as f:
        return json.load(f)


_REGISTRY = load_registry()


def normalize_assay_field_format(df: pd.DataFrame) -> None:
    """Respell each Assay free-text column to its registry spelling where a
    cell differs from it only in case or spacing (auto_fix, to correct in
    the Excel file too)."""
    for column in ASSAY_FORMAT_COLUMNS:
        by_key = _REGISTRY.get(column)
        if not by_key:
            continue
        for idx, raw in df[ASSAY][column].items():
            if not isinstance(raw, str) or not raw.strip():
                continue
            canonical = by_key.get(registry_key(raw))
            if not canonical or canonical == raw:
                continue
            emit(
                Issue(
                    severity="auto_fix",
                    code="assay_field_format",
                    message="Value differs from its registry spelling only in "
                    "case or spacing; the export writes the registry spelling. "
                    "Correct it in the Excel file.",
                    group=ASSAY,
                    column=column,
                    rows=[idx],
                    value=raw,
                    suggested_value=canonical,
                )
            )
            df.loc[idx, (ASSAY, column)] = canonical
