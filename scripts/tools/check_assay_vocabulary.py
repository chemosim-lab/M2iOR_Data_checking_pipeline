# pipeline/scripts/tools/check_assay_vocabulary.py
"""List the Assay values of Excel files that are outside the allowed lists.

Reads each study the way the pipeline does (strip, lowercase, aliases,
respelling to the listed spelling) and reports, per column, the values that
are not in `assay_vocabulary.json` (they block the export) and the pending
ones (they only warn), with the studies and row counts. Offline and
read-only: use it to see what a new delivery or a corpus holds before a
value is added to a list.

Usage:
    uv run python -m scripts.tools.check_assay_vocabulary                # every study under input/
    uv run python -m scripts.tools.check_assay_vocabulary input/Pelz_2006.xlsx
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from scripts.normalize_dataframe import normalize_df
from scripts.process_assay import (
    _VOCABULARY,
    normalize_assay_field_format,
    registry_key,
    unlisted_values,
)
from scripts.read_excel import get_raw_data_from_excel_file

INPUT_DIR = Path(__file__).resolve().parent.parent.parent / "input"


def main(argv: list[str]) -> int:
    paths = [Path(a) for a in argv] or sorted(INPUT_DIR.glob("*.xlsx"))
    report: dict[str, dict[str, dict[str, dict[str, int]]]] = {}
    blocking = 0
    for path in paths:
        df = normalize_df(get_raw_data_from_excel_file(path))
        normalize_assay_field_format(df)
        for column, values in unlisted_values(df).items():
            for value, rows in values.items():
                entry = _VOCABULARY[column].entries.get(registry_key(value))
                status = "pending" if entry is not None else "not_allowed"
                blocking += entry is None
                report.setdefault(status, {}).setdefault(column, {}).setdefault(
                    value, {}
                )[path.stem] = len(rows)
    json.dump(report, sys.stdout, ensure_ascii=False, indent=2)
    print()
    return 1 if blocking else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
