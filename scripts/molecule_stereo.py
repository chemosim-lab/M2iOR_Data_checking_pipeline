# pipeline/scripts/molecule_stereo.py
"""Stereochemistry ambiguity classification and structure rendering.

Classifies a molecule's SMILES as:
  - "monomolecular"  -> stereochemistry is fully defined (or none exists)
  - "sum of isomers" -> at least one stereo element (R/S or Z/E) is undefined

and renders a PNG of the structure with undefined stereocenters highlighted
in red and undefined double bonds highlighted in blue.
"""

import io
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

MONOMOLECULAR = "monomolecular"
SUM_OF_ISOMERS = "sum of isomers"

_UNDEFINED_ATOM_COLOR = (1.0, 0.4, 0.4)
_UNDEFINED_BOND_COLOR = (0.4, 0.6, 1.0)
_IMAGE_SIZE = (300, 260)

_LEGEND_FONT_SIZE = 13
_LEGEND_LINE_HEIGHT = 18
_LEGEND_PADDING = 6
_LEGEND_TEXT_COLOR = (70, 70, 70)


def _to_rgb255(color: tuple[float, float, float]) -> tuple[int, int, int]:
    r, g, b = color
    return (round(r * 255), round(g * 255), round(b * 255))

# Modern (non-legacy) stereo perception correctly handles ring-size
# constraints, terminal double bonds, and symmetric substituents (e.g. two
# identical methyl groups on a double bond are NOT stereogenic).
Chem.SetUseLegacyStereoPerception(False)  # noqa: FBT003 - third-party API

# SANITIZE_SETAROMATICITY is mandatory here: without it, aromatic rings
# (e.g. c1ccccc1 or C1=CC=CC=C1) are kept as alternating localized double
# bonds and RDKit would wrongly flag them as undefined E/Z double bonds.
_SANITIZE_FLAGS = (
    Chem.SANITIZE_CLEANUP
    | Chem.SANITIZE_PROPERTIES
    | Chem.SANITIZE_SYMMRINGS
    | Chem.SANITIZE_KEKULIZE
    | Chem.SANITIZE_SETAROMATICITY
    | Chem.SANITIZE_SETCONJUGATION
    | Chem.SANITIZE_SETHYBRIDIZATION
    | Chem.SANITIZE_CLEANUPCHIRALITY
    | Chem.SANITIZE_ADJUSTHS
)


@dataclass
class StereoResult:
    status: str | None  # None means the SMILES could not be parsed
    mol: Chem.Mol | None
    undefined_atoms: list[int]
    undefined_bonds: list[int]


def analyze_stereo(smiles: str) -> StereoResult:
    """Classify a SMILES and report its undefined stereo elements."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return StereoResult(None, None, [], [])

    try:
        Chem.SanitizeMol(mol, sanitizeOps=_SANITIZE_FLAGS)
    except Exception:  # noqa: BLE001 - any sanitization failure means unusable mol
        return StereoResult(None, None, [], [])

    Chem.AssignStereochemistry(
        mol, cleanIt=True, force=True, flagPossibleStereoCenters=True
    )

    undefined_atoms: list[int] = []
    undefined_bonds: list[int] = []
    for element in Chem.FindPotentialStereo(mol):
        if element.specified != Chem.StereoSpecified.Unspecified:
            continue
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            undefined_atoms.append(element.centeredOn)
        elif element.type == Chem.StereoType.Bond_Double:
            undefined_bonds.append(element.centeredOn)

    status = SUM_OF_ISOMERS if (undefined_atoms or undefined_bonds) else MONOMOLECULAR
    return StereoResult(status, mol, undefined_atoms, undefined_bonds)


def _draw_legend_line(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    dot_color: tuple[int, int, int],
    label: str,
    rest: str,
    font: ImageFont.ImageFont,
) -> None:
    x, y = xy
    dot_radius = _LEGEND_FONT_SIZE // 3
    dot_cy = y + _LEGEND_FONT_SIZE // 2
    draw.ellipse(
        (x, dot_cy - dot_radius, x + 2 * dot_radius, dot_cy + dot_radius),
        fill=dot_color,
    )
    x += 2 * dot_radius + 4
    draw.text((x, y), label, font=font, fill=dot_color)
    x += round(draw.textlength(label, font=font))
    draw.text((x, y), rest, font=font, fill=_LEGEND_TEXT_COLOR)


def _append_legend(
    png_bytes: bytes, legend_lines: list[tuple[tuple[int, int, int], str, str]]
) -> bytes:
    """Append a white strip below the image explaining the highlight colors."""
    base = Image.open(io.BytesIO(png_bytes)).convert("RGB")
    font = ImageFont.load_default(size=_LEGEND_FONT_SIZE)

    strip_height = _LEGEND_PADDING * 2 + _LEGEND_LINE_HEIGHT * len(legend_lines)
    canvas = Image.new("RGB", (base.width, base.height + strip_height), "white")
    canvas.paste(base, (0, 0))

    draw = ImageDraw.Draw(canvas)
    y = base.height + _LEGEND_PADDING
    for dot_color, label, rest in legend_lines:
        _draw_legend_line(draw, (_LEGEND_PADDING, y), dot_color, label, rest, font)
        y += _LEGEND_LINE_HEIGHT

    buf = io.BytesIO()
    canvas.save(buf, format="PNG")
    return buf.getvalue()


def render_stereo_image(result: StereoResult) -> bytes:
    """Render the molecule as a PNG, highlighting undefined stereo elements.

    The status label is drawn under the structure in sentence case
    (e.g. "Monomolecular", "Sum of isomers"). If any stereo element is
    undefined, a legend explaining the highlight colors is appended below.
    """
    if result.mol is None or result.status is None:
        msg = "Cannot render an image for an unparseable SMILES"
        raise ValueError(msg)

    atom_colors = dict.fromkeys(result.undefined_atoms, _UNDEFINED_ATOM_COLOR)
    bond_colors = dict.fromkeys(result.undefined_bonds, _UNDEFINED_BOND_COLOR)

    drawer = rdMolDraw2D.MolDraw2DCairo(*_IMAGE_SIZE)
    drawer.drawOptions().addStereoAnnotation = True
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        result.mol,
        legend=result.status.capitalize(),
        highlightAtoms=result.undefined_atoms,
        highlightBonds=result.undefined_bonds,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()

    legend_lines: list[tuple[tuple[int, int, int], str, str]] = []
    if result.undefined_atoms:
        legend_lines.append(
            (_to_rgb255(_UNDEFINED_ATOM_COLOR), "red atom", " = undefined R/S center")
        )
    if result.undefined_bonds:
        legend_lines.append(
            (_to_rgb255(_UNDEFINED_BOND_COLOR), "blue bond", " = undefined Z/E bond")
        )

    if not legend_lines:
        return png_bytes
    return _append_legend(png_bytes, legend_lines)


def export_stereo_classification(
    unique_cids: list[int],
    smiles_map: dict[int, str],
    images_dir: Path,
    *,
    overwrite: bool = False,
) -> dict[int, str]:
    """Classify each CID's SMILES and save its highlighted structure image.

    Images are saved as "<cid>.png" in *images_dir* so the front-end can
    associate a CID with its rendered structure. A CID whose image already
    exists is not re-rendered unless *overwrite* is set. Returns a cid ->
    status map for CIDs whose SMILES could be parsed.
    """
    images_dir.mkdir(parents=True, exist_ok=True)
    mixture_map: dict[int, str] = {}

    for cid in unique_cids:
        smiles = smiles_map.get(cid)
        if not smiles:
            continue

        result = analyze_stereo(smiles)
        if result.status is None:
            continue
        mixture_map[cid] = result.status

        image_path = images_dir / f"{cid}.png"
        if overwrite or not image_path.exists():
            image_path.write_bytes(render_stereo_image(result))

    return mixture_map
