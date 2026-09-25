# pipeline/scripts/report.py
"""Structured validation issues and the machine-readable check report.

Validation steps raise `ValidationError` - a `ValueError`, so existing
`except ValueError` handling keeps working - carrying one `Issue` per
problem, and call `emit()` for non-blocking findings (warnings, values the
pipeline silently rewrites, informational enrichments). `IssueCollector`
runs the pipeline stage by stage and gathers both, and `build_report` turns
them into the document printed by `main.py --check --report json`.

Only "error" issues block the export; the other severities never change
what the pipeline accepts, they only surface what it does.
"""

from __future__ import annotations

import json
import math
import traceback
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Literal

import pandas as pd
from openpyxl.utils import get_column_letter

from scripts.columns import ACCESSION, GENE_NAME, RECEPTOR_NAME, UNIPROT_ID

REPORT_SCHEMA_VERSION = 1

# Two header rows (group labels, column names) precede the data, so the row
# at DataFrame index 0 is Excel row 3 (see scripts/read_excel.py).
EXCEL_HEADER_ROWS = 2
EXCEL_FIRST_DATA_ROW = EXCEL_HEADER_ROWS + 1

Severity = Literal["error", "warning", "auto_fix", "info"]
SEVERITY_ORDER: tuple[Severity, ...] = ("error", "warning", "auto_fix", "info")

StageStatus = Literal["passed", "failed", "crashed", "skipped"]
_STAGE_STATUS_RANK: dict[StageStatus, int] = {
    "passed": 0,
    "skipped": 0,
    "failed": 1,
    "crashed": 2,
}

# Columns the pipeline renames after reading (see main.py), mapped back to the
# header the Excel file actually uses, so issues point at the right header.
_EXCEL_HEADER_BY_PIPELINE_COLUMN = {
    RECEPTOR_NAME: GENE_NAME,
    ACCESSION: UNIPROT_ID,
}

# Beyond this many distinct raw values, `excel_values` is truncated.
_MAX_EXCEL_VALUES = 20


def excel_row(index: int) -> int:
    """1-based Excel row of a DataFrame row index on the data sheet."""
    return int(index) + EXCEL_FIRST_DATA_ROW


@dataclass
class Issue:
    """One finding about the input Excel file.

    `rows` are DataFrame row indices (converted to Excel rows in the report).
    `cells` is only set directly for findings outside the data rows, such as
    a header cell. Messages must not mention row numbers: issues that only
    differ by their rows are merged into one in the report.
    """

    code: str
    message: str
    severity: Severity = "error"
    group: str | None = None
    column: str | None = None
    rows: list[int] = field(default_factory=list[int])
    value: Any = None
    suggested_value: Any = None
    details: dict[str, Any] = field(default_factory=dict[str, Any])
    cells: list[str] | None = None
    stage: str | None = None


class ValidationError(ValueError):
    """A blocking validation failure, carrying its structured issues.

    `str(error)` stays the human-readable message the pipeline has always
    printed and stored in registry.json."""

    def __init__(self, message: str, issues: list[Issue]) -> None:
        super().__init__(message)
        self.issues = issues


_active_collector: ContextVar[IssueCollector | None] = ContextVar(
    "active_issue_collector", default=None
)


def emit(issue: Issue) -> None:
    """Record a non-blocking finding with the active collector, if any.

    A no-op outside `IssueCollector.stage()`, so library code can call it
    unconditionally."""
    collector = _active_collector.get()
    if collector is not None:
        collector.add(issue)


class IssueCollector:
    """Gathers issues stage by stage instead of stopping at the first error.

    An exception escaping a `stage()` block aborts the rest of that stage
    only; `check()` runs one independent step and carries on either way.
    Unexpected (non-`ValueError`) exceptions are recorded as an
    `internal_error` when `catch_unexpected` is set, re-raised otherwise.
    """

    def __init__(self, *, catch_unexpected: bool) -> None:
        self.issues: list[Issue] = []
        self.stages: dict[str, StageStatus] = {}
        # str() of every blocking failure, in order - the registry message.
        self.failure_messages: list[str] = []
        self._catch_unexpected = catch_unexpected
        self._stage: str | None = None

    @property
    def has_errors(self) -> bool:
        return any(issue.severity == "error" for issue in self.issues)

    @property
    def crashed(self) -> bool:
        return "crashed" in self.stages.values()

    def add(self, issue: Issue) -> None:
        if issue.stage is None:
            issue.stage = self._stage
        self.issues.append(issue)
        if issue.severity == "error" and self._stage is not None:
            self._mark(self._stage, "failed")

    def skip(self, stage: str) -> None:
        self.stages.setdefault(stage, "skipped")

    @contextmanager
    def stage(self, name: str) -> Iterator[None]:
        previous = self._stage
        self._stage = name
        self.stages.setdefault(name, "passed")
        token = _active_collector.set(self)
        try:
            yield
        except Exception as e:  # noqa: BLE001 - recorded, see _record_failure
            self._record_failure(e)
        finally:
            _active_collector.reset(token)
            self._stage = previous

    def check(self, step: Callable[..., object], *args: Any, **kwargs: Any) -> bool:
        """Run one independent validation step; return whether it passed."""
        try:
            step(*args, **kwargs)
        except Exception as e:  # noqa: BLE001 - recorded, see _record_failure
            self._record_failure(e)
            return False
        return True

    def _record_failure(self, error: Exception) -> None:
        if isinstance(error, ValidationError):
            self.failure_messages.append(str(error))
            for issue in error.issues:
                issue.severity = "error"
                self.add(issue)
        elif isinstance(error, ValueError):
            self.failure_messages.append(str(error))
            self.add(Issue(code="unstructured_error", message=str(error)))
        else:
            if not self._catch_unexpected:
                raise error
            message = f"{type(error).__name__}: {error}"
            self.failure_messages.append(f"internal error: {message}")
            self.add(
                Issue(
                    code="internal_error",
                    message=message,
                    details={"traceback": traceback.format_exc()},
                )
            )
            if self._stage is not None:
                self._mark(self._stage, "crashed")

    def _mark(self, stage: str, status: StageStatus) -> None:
        current = self.stages.get(stage, "passed")
        if _STAGE_STATUS_RANK[status] > _STAGE_STATUS_RANK[current]:
            self.stages[stage] = status


@dataclass
class ExcelLayout:
    """Where the data sheet's columns sit, and its values as first read.

    `raw` is the DataFrame exactly as read from the Excel file, before any
    normalization, sharing the processed DataFrame's row index."""

    sheet: str
    raw: pd.DataFrame
    _positions: dict[tuple[str, str], int] = field(
        default_factory=dict[tuple[str, str], int]
    )

    @classmethod
    def from_raw(cls, raw: pd.DataFrame) -> ExcelLayout:
        positions: dict[tuple[str, str], int] = {}
        for pos, (group, column) in enumerate(raw.columns):
            positions.setdefault((str(group), str(column)), pos)
        return cls(
            sheet=str(raw.attrs.get("sheet_name", "")),
            raw=raw.copy(deep=True),
            _positions=positions,
        )

    @staticmethod
    def excel_header(column: str) -> str:
        return _EXCEL_HEADER_BY_PIPELINE_COLUMN.get(column, column)

    def position(self, group: str | None, column: str | None) -> int | None:
        if group is None or column is None:
            return None
        return self._positions.get((group, self.excel_header(column)))

    def excel_column(self, group: str | None, column: str | None) -> str | None:
        pos = self.position(group, column)
        return get_column_letter(pos + 1) if pos is not None else None

    def raw_values(self, group: str | None, column: str | None, rows: list[int]) -> list[Any]:
        """Distinct values of the given rows' cells, as read from the Excel file."""
        pos = self.position(group, column)
        if pos is None or not rows:
            return []
        series = self.raw.iloc[:, pos]
        values: list[Any] = []
        for row in rows:
            if row not in series.index:
                continue
            value = _jsonable(series.loc[row])
            if value not in values:
                values.append(value)
            if len(values) >= _MAX_EXCEL_VALUES:
                break
        return values


def _jsonable(value: Any) -> Any:
    """Convert pandas/numpy scalars and containers into plain JSON values."""
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}  # pyright: ignore[reportUnknownVariableType]
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]  # pyright: ignore[reportUnknownVariableType]
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return None if math.isnan(value) else value
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):  # numpy scalar
        return _jsonable(value.item())
    return str(value)


def _merge_key(issue: Issue) -> str:
    return json.dumps(
        [
            issue.severity,
            issue.stage,
            issue.code,
            issue.group,
            issue.column,
            issue.message,
            _jsonable(issue.value),
            _jsonable(issue.suggested_value),
            _jsonable(issue.details),
            issue.cells,
        ],
        sort_keys=True,
        ensure_ascii=False,
    )


def merge_issues(issues: list[Issue]) -> list[Issue]:
    """Merge issues differing only by their rows (e.g. one wrong molecule
    name repeated on every row testing it), keeping first-seen order."""
    merged: dict[str, Issue] = {}
    for issue in issues:
        key = _merge_key(issue)
        if key in merged:
            merged[key].rows.extend(issue.rows)
        else:
            merged[key] = Issue(
                code=issue.code,
                message=issue.message,
                severity=issue.severity,
                group=issue.group,
                column=issue.column,
                rows=list(issue.rows),
                value=issue.value,
                suggested_value=issue.suggested_value,
                details=issue.details,
                cells=issue.cells,
                stage=issue.stage,
            )
    for issue in merged.values():
        issue.rows = sorted(set(issue.rows))
    return list(merged.values())


def _serialize_issue(issue: Issue, layout: ExcelLayout | None) -> dict[str, Any]:
    letter = layout.excel_column(issue.group, issue.column) if layout else None
    rows = [excel_row(r) for r in issue.rows]
    cells = issue.cells
    if cells is None and letter is not None:
        cells = [f"{letter}{row}" for row in rows]
    return {
        "severity": issue.severity,
        "stage": issue.stage,
        "code": issue.code,
        "message": issue.message,
        "sheet": layout.sheet if layout else None,
        "group": issue.group,
        "column": ExcelLayout.excel_header(issue.column) if issue.column else None,
        "excel_column": letter,
        "excel_rows": rows,
        "cells": cells or [],
        "value": _jsonable(issue.value),
        "excel_values": layout.raw_values(issue.group, issue.column, issue.rows)
        if layout
        else [],
        "suggested_value": _jsonable(issue.suggested_value),
        "details": _jsonable(issue.details),
    }


def build_report(
    collector: IssueCollector,
    *,
    mode: Literal["check", "export"],
    study_id: str,
    excel_path: str,
    file_sha256: str,
    layout: ExcelLayout | None,
    blast_mode: str,
    exported_csv: str | None = None,
) -> dict[str, Any]:
    stage_order = list(collector.stages)
    issues = merge_issues(collector.issues)
    issues.sort(
        key=lambda i: (
            SEVERITY_ORDER.index(i.severity),
            stage_order.index(i.stage) if i.stage in stage_order else len(stage_order),
            i.rows[0] if i.rows else -1,
        )
    )
    summary = {severity: 0 for severity in SEVERITY_ORDER}
    for issue in issues:
        summary[issue.severity] += 1

    if collector.crashed:
        status = "crashed"
    elif collector.has_errors:
        status = "failed"
    else:
        status = "passed"

    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "mode": mode,
        "study_id": study_id,
        "excel_path": excel_path,
        "file_sha256": file_sha256,
        "sheet": layout.sheet if layout else None,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "blast_mode": blast_mode,
        "status": status,
        "stages": dict(collector.stages),
        "summary": summary,
        "exported_csv": exported_csv,
        "issues": [_serialize_issue(issue, layout) for issue in issues],
    }


def _format_location(issue: dict[str, Any]) -> str:
    items = issue["cells"] or [str(row) for row in issue["excel_rows"]]
    label = "cells" if issue["cells"] else "rows"
    if not items:
        return f"{label}: -"
    if len(items) <= 4:
        return f"{label}: {','.join(items)}"
    return f"{label}: {','.join(items[:3])} (+{len(items) - 3} more)"


def format_text_report(report: dict[str, Any]) -> str:
    """Human-readable rendering of `build_report`'s output. "info" issues
    are only counted per code; the JSON report lists them all."""
    counts = ", ".join(f"{n} {sev}" for sev, n in report["summary"].items() if n)
    lines = [
        f"{report['study_id']} [{report['mode']}] - {report['status'].upper()}"
        f" ({counts or 'no issues'}) - sheet: {report['sheet']}",
        "stages: "
        + ", ".join(f"{name}={status}" for name, status in report["stages"].items()),
    ]
    if report["exported_csv"]:
        lines.append(f"exported: {report['exported_csv']}")
    info_counts: dict[str, int] = {}
    for issue in report["issues"]:
        if issue["severity"] == "info":
            info_counts[issue["code"]] = info_counts.get(issue["code"], 0) + 1
            continue
        where = "/".join(p for p in (issue["group"], issue["column"]) if p)
        lines.append(
            f"[{issue['severity'].upper()}] {issue['stage']}/{issue['code']}"
            f"  {where or '-'}  {_format_location(issue)}"
            + (f"  value: {issue['value']!r}" if issue["value"] is not None else "")
        )
        lines.append(f"    {issue['message']}")
        if issue["suggested_value"] is not None:
            lines.append(f"    suggested: {issue['suggested_value']!r}")
    for code, count in info_counts.items():
        lines.append(f"[INFO] {code}: {count} issue(s), see the JSON report")
    return "\n".join(lines)
