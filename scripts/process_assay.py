# pipeline/scripts/process_assay.py
"""Allowed values for the free-text Assay columns (`assay_vocabulary.json`).

Each column has a closed list of allowed values, a definition, and the
values that are still in the curated files but not validated (`pending`):

- a cell that differs from a listed value only in case, spacing or the
  delta sign is respelled to the listed spelling (`assay_field_format`,
  auto_fix, to correct in the Excel file too);
- a pending value is kept and reported (`assay_value_pending`, warning);
- any other value blocks the export (`value_not_allowed`) with the closest
  listed value as a suggestion: a new instrument, vector or neuron class is
  added to the list on purpose, with its definition, not by accident;
- a listed value that declares what it is used with (`requires`) must sit
  on a row whose technique or expression system is one of those
  (`assay_inconsistent`): an amplifier on a fluorescence row, or a plasmid
  as the driver of a Drosophila neuron, blocks the export.

Experimental technique and the unit columns keep their own lists in
`scripts/process_responses.py`.
"""

import difflib
import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from scripts.columns import ASSAY
from scripts.report import Issue, ValidationError, emit

VOCABULARY_PATH = Path(__file__).resolve().parent.parent / "assay_vocabulary.json"

# Greek small delta and the increment sign are folded together, so that
# "Δ fluorescence" typed with the Greek letter matches "∆ fluorescence".
_DELTA_FOLD = {0x3B4: 0x2206}


def registry_key(value: str) -> str:
    """Normalized form used to match a cell against the vocabulary: case,
    surrounding/repeated whitespace and the delta sign only, so two values
    differing in substance never collide."""
    return " ".join(value.strip().lower().translate(_DELTA_FOLD).split())


@dataclass(frozen=True)
class Entry:
    canonical: str
    pending: bool = False
    note: str | None = None
    requires: dict[str, list[str]] = field(default_factory=dict[str, list[str]])


@dataclass(frozen=True)
class ColumnVocabulary:
    definition: str
    entries: dict[str, Entry]  # registry key -> entry

    @property
    def allowed(self) -> list[str]:
        return sorted(e.canonical for e in self.entries.values() if not e.pending)


Vocabulary = dict[str, ColumnVocabulary]


def parse_vocabulary(raw: dict[str, dict[str, object]]) -> Vocabulary:
    vocabulary: Vocabulary = {}
    for column, spec in raw.items():
        entries: dict[str, Entry] = {}
        allowed: dict[str, dict[str, object]] = spec.get("allowed", {})  # type: ignore[assignment]
        pending: dict[str, str] = spec.get("pending", {})  # type: ignore[assignment]
        for value, info in allowed.items():
            entries[registry_key(value)] = Entry(
                canonical=value,
                note=info.get("note"),  # type: ignore[arg-type]
                requires=info.get("requires", {}),  # type: ignore[arg-type]
            )
        for value, note in pending.items():
            entries.setdefault(
                registry_key(value), Entry(canonical=value, pending=True, note=note)
            )
        vocabulary[column] = ColumnVocabulary(
            definition=str(spec.get("definition", "")), entries=entries
        )
    return vocabulary


def load_vocabulary(path: Path = VOCABULARY_PATH) -> Vocabulary:
    with path.open(encoding="utf-8") as f:
        return parse_vocabulary(json.load(f))


_VOCABULARY = load_vocabulary()


def _cell_groups(df: pd.DataFrame, column: str) -> dict[str, list[int]]:
    """Non-empty cells of an Assay column, grouped by their (stripped) text."""
    groups: dict[str, list[int]] = {}
    for idx, raw in df[ASSAY][column].items():
        if pd.isna(raw):
            continue
        text = str(raw).strip()
        if text:
            groups.setdefault(text, []).append(idx)
    return groups


def normalize_assay_field_format(df: pd.DataFrame) -> None:
    """Respell each cell to its listed spelling where it differs only in case,
    spacing or delta sign (auto_fix), and report the pending values."""
    for column, vocabulary in _VOCABULARY.items():
        if (ASSAY, column) not in df.columns:
            continue
        for raw, rows in _cell_groups(df, column).items():
            entry = vocabulary.entries.get(registry_key(raw))
            if entry is None:
                continue
            if entry.pending:
                emit(
                    Issue(
                        severity="warning",
                        code="assay_value_pending",
                        message=f"'{raw}' is not validated: {entry.note}",
                        group=ASSAY,
                        column=column,
                        rows=rows,
                        value=raw,
                    )
                )
            if entry.canonical != raw:
                emit(
                    Issue(
                        severity="auto_fix",
                        code="assay_field_format",
                        message="Value differs from its listed spelling only in "
                        "case, spacing or delta sign; the export writes the "
                        "listed spelling. Correct it in the Excel file.",
                        group=ASSAY,
                        column=column,
                        rows=rows,
                        value=raw,
                        suggested_value=entry.canonical,
                    )
                )
                for idx in rows:
                    df.loc[idx, (ASSAY, column)] = entry.canonical


def validate_assay_vocabulary(df: pd.DataFrame) -> None:
    """Raise ValidationError for every value outside the allowed lists."""
    issues: list[Issue] = []
    for column, vocabulary in _VOCABULARY.items():
        if (ASSAY, column) not in df.columns:
            continue
        allowed = vocabulary.allowed
        by_key = {registry_key(v): v for v in allowed}
        for raw, rows in _cell_groups(df, column).items():
            if registry_key(raw) in vocabulary.entries:
                continue
            closest = difflib.get_close_matches(
                registry_key(raw), list(by_key), n=1, cutoff=0.6
            )
            issues.append(
                Issue(
                    code="value_not_allowed",
                    message=f"'{raw}' is not an allowed '{column}' value.",
                    group=ASSAY,
                    column=column,
                    rows=rows,
                    value=raw,
                    suggested_value=by_key[closest[0]] if closest else None,
                    details={"allowed": allowed, "definition": vocabulary.definition},
                )
            )
    if issues:
        bad = sorted({f"{i.column}: {i.value!r}" for i in issues})
        raise ValidationError(f"Values outside the allowed lists: {bad}", issues)


def validate_assay_consistency(df: pd.DataFrame) -> None:
    """Raise ValidationError when a value sits on a row whose technique or
    expression system it is not used with (`requires` in the vocabulary)."""
    issues: list[Issue] = []
    for column, vocabulary in _VOCABULARY.items():
        if (ASSAY, column) not in df.columns:
            continue
        for raw, rows in _cell_groups(df, column).items():
            entry = vocabulary.entries.get(registry_key(raw))
            if entry is None or entry.pending:
                continue
            for other, expected in entry.requires.items():
                if (ASSAY, other) not in df.columns:
                    continue
                expected_keys = {registry_key(v) for v in expected}
                by_other: dict[str, list[int]] = {}
                for idx in rows:
                    cell = df[ASSAY][other].loc[idx]
                    if pd.isna(cell) or not str(cell).strip():
                        continue
                    if registry_key(str(cell)) not in expected_keys:
                        by_other.setdefault(str(cell).strip(), []).append(idx)
                for other_value, other_rows in by_other.items():
                    issues.append(
                        Issue(
                            code="assay_inconsistent",
                            message=f"'{raw}' is not used with {other} "
                            f"'{other_value}' (expected: {', '.join(expected)}).",
                            group=ASSAY,
                            column=column,
                            rows=other_rows,
                            value=raw,
                            details={"column": other, "found": other_value, "expected": expected},
                        )
                    )
    if issues:
        bad = sorted({f"{i.column}: {i.value!r}" for i in issues})
        raise ValidationError(f"Values inconsistent with the rest of the row: {bad}", issues)


def unlisted_values(df: pd.DataFrame) -> dict[str, dict[str, list[int]]]:
    """Column -> value -> rows, for the values outside the allowed lists,
    pending ones included. Used by `scripts/tools/check_assay_vocabulary.py`."""
    found: dict[str, dict[str, list[int]]] = {}
    for column, vocabulary in _VOCABULARY.items():
        if (ASSAY, column) not in df.columns:
            continue
        for raw, rows in _cell_groups(df, column).items():
            entry = vocabulary.entries.get(registry_key(raw))
            if entry is None or entry.pending:
                found.setdefault(column, {})[raw] = rows
    return found
