# pipeline/scripts/tools/build_assay_registry.py
"""(Re)build assay_registry.json from the current, curated Excel corpus.

Scans every study under input/ for the Assay free-text columns (see
scripts/process_assay.py), collects their values as the pipeline would see
them (after normalize_df's strip/lowercase), and writes one canonical
spelling per registry key (case/spacing-insensitive). Run this again after a
legitimate new equipment/construct name has been added to a curated study,
so the registry keeps learning; it is not run automatically by the pipeline.

Two Excel values sharing a registry key but spelled differently would be a
real inconsistency in the "cleaned" corpus: reported and kept out of the
registry rather than picking one arbitrarily.

Usage:
    uv run python -m scripts.tools.build_assay_registry
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts.columns import ASSAY
from scripts.normalize_dataframe import normalize_df
from scripts.process_assay import ASSAY_FORMAT_COLUMNS, REGISTRY_PATH, registry_key
from scripts.read_excel import get_raw_data_from_excel_file

INPUT_DIR = Path(__file__).resolve().parent.parent.parent / "input"


def build_registry(input_dir: Path = INPUT_DIR) -> dict[str, dict[str, str]]:
    registry: dict[str, dict[str, str]] = {column: {} for column in ASSAY_FORMAT_COLUMNS}
    conflicting_keys: set[tuple[str, str]] = set()

    for path in sorted(input_dir.glob("*.xlsx")):
        df = normalize_df(get_raw_data_from_excel_file(path))
        for column in ASSAY_FORMAT_COLUMNS:
            if (ASSAY, column) not in df.columns:
                continue
            for raw in df[ASSAY][column].dropna():
                value = str(raw).strip()
                if not value:
                    continue
                key = registry_key(value)
                existing = registry[column].get(key)
                if existing is None:
                    registry[column][key] = value
                elif existing != value:
                    conflicting_keys.add((column, key))
                    print(
                        f"conflict, kept out of the registry: {path.name}: "
                        f"'{column}' has both {existing!r} and {value!r} for "
                        f"registry key {key!r}"
                    )

    for column, key in conflicting_keys:
        registry[column].pop(key, None)

    return {column: dict(sorted(values.items())) for column, values in registry.items()}


def main() -> None:
    registry = build_registry()
    with REGISTRY_PATH.open("w", encoding="utf-8") as f:
        json.dump(registry, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write("\n")
    total = sum(len(values) for values in registry.values())
    print(f"wrote {REGISTRY_PATH} ({total} entries across {len(registry)} columns)")


if __name__ == "__main__":
    main()
