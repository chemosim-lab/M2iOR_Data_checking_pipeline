# pipeline/scripts/tools/backfill_stereo_images.py
"""Backfill per-CID stereo-ambiguity structure images from the local PubChem cache.

Scans every cached PubChem molecule record (./cache/molecules) and, for each
CID with a usable SMILES, classifies its stereochemistry ambiguity and writes
a highlighted structure PNG named "<cid>.png" to the target images directory.
CIDs whose image already exists are left untouched.

This does not touch already-exported CSVs or their Mixture column -- rerun
the main pipeline (main.py --force) to refresh those.

Usage:
    python -m scripts.tools.backfill_stereo_images
    python -m scripts.tools.backfill_stereo_images --images-dir path/to/images
"""

import argparse
import json
from pathlib import Path
from typing import Any

from scripts.molecule_stereo import export_stereo_classification

_CACHE_DIR = Path("./cache/molecules")
# The M2iOR website project is expected next to this pipeline's directory
_DEFAULT_IMAGES_DIR = (
    Path(__file__).resolve().parents[3]
    / "M2iOR_web_public-main"
    / "public"
    / "images"
    / "molecules"
)


def _extract_smiles(props: list[dict[str, Any]]) -> str | None:
    isomeric: str | None = None
    absolute: str | None = None
    for prop in props:
        urn = prop.get("urn", {})
        if urn.get("label") != "SMILES":
            continue
        if urn.get("name") == "Isomeric":
            isomeric = prop.get("value", {}).get("sval")
        elif urn.get("name") == "Absolute":
            absolute = prop.get("value", {}).get("sval")
    return isomeric or absolute


def load_cached_smiles(cache_dir: Path) -> dict[int, str]:
    """Read every cached PubChem record and return a cid -> SMILES map."""
    smiles_map: dict[int, str] = {}
    for path in cache_dir.glob("*.json"):
        data = json.loads(path.read_text(encoding="utf-8"))
        compounds = data.get("PC_Compounds", [])
        if not compounds:
            continue
        cid = compounds[0].get("id", {}).get("id", {}).get("cid")
        smiles = _extract_smiles(compounds[0].get("props", []))
        if cid and smiles:
            smiles_map[cid] = smiles
    return smiles_map


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Backfill stereo-ambiguity structure images from the local PubChem cache."
        )
    )
    _ = parser.add_argument("--cache-dir", type=Path, default=_CACHE_DIR)
    _ = parser.add_argument("--images-dir", type=Path, default=_DEFAULT_IMAGES_DIR)
    _ = parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-render images that already exist (e.g. after a rendering change)",
    )
    args = parser.parse_args()

    smiles_map = load_cached_smiles(args.cache_dir)
    print(f"Loaded {len(smiles_map)} cached CIDs with a usable SMILES")

    mixture_map = export_stereo_classification(
        list(smiles_map), smiles_map, args.images_dir, overwrite=args.overwrite
    )
    print(f"Classified {len(mixture_map)} CIDs, images saved to {args.images_dir}")


if __name__ == "__main__":
    main()
