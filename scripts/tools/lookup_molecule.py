# pipeline/scripts/tools/lookup_molecule.py
"""Look up a molecule and diagnose an Excel row's Name/CID/CAS consistency.

Uses the pipeline's own sources, cache and rules (scripts/fetch_data.py,
scripts/process_molecules.py): PubChem records and synonyms, CAS Common
Chemistry, CAS -> CID resolution, and the same CID/CAS reconciliation and
Molecule Name check as main.py. For a row it reports:

  - "pipeline": how main.py treats the row as is (errors, silent rewrites),
  - "compounds" / "cas": what each CID and CAS number found along the way is,
  - "lookups": which CIDs PubChem finds for the name / SMILES / InChIKey,
  - "suggestions": candidate corrections of the Excel row, each re-checked
    with the pipeline's rules ("pipeline_after") - candidates, not decisions.

Run from the pipeline root (the cache lives in ./cache). Prints JSON.

Usage:
    uv run python -m scripts.tools.lookup_molecule --name "4-Methylcyclohexanol" --cid 87468080 --cas 25639-42-3
    uv run python -m scripts.tools.lookup_molecule --smiles "CC1CCC(CC1)O"
    uv run python -m scripts.tools.lookup_molecule --excel input/Pelz_2006.xlsx --row 8 --row 17
    uv run python -m scripts.tools.lookup_molecule --report check.json
"""

from __future__ import annotations

import argparse
import difflib
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.columns import CAS, CID, INCHIKEY, MIXTURE, MOLECULE, MOLECULE_NAME, SMILES
from scripts.fetch_data import (
    fetch_cas_common_chemistry_details,
    fetch_cas_to_cid_map,
    fetch_cids_by_inchikey,
    fetch_cids_by_name,
    fetch_cids_by_smiles,
    fetch_pubchem_data,
)
from scripts.molecule_stereo import analyze_stereo
from scripts.process_molecules import (
    _normalize_name,
    canonical_name,
    registry_spelling,
    _parse_cid_cache,
    check_cid_cas,
    molecule_name_matches,
    normalize_dashes,
    parse_cas_cell,
    parse_cid_cell,
)
from scripts.read_excel import get_raw_data_from_excel_file
from scripts.report import (
    EXCEL_FIRST_DATA_ROW,
    ExcelLayout,
    IssueCollector,
    _jsonable,
)

_PUBCHEM_COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
_CAS_COMMON_CHEMISTRY_DETAIL_URL = "https://commonchemistry.cas.org/detail?cas_rn={cas}"

# Check report issues this tool can diagnose (see --report).
MOLECULE_ISSUE_CODES = {
    "cid_and_cas_missing",
    "cid_cas_name_conflict",
    "molecule_name_mismatch",
    "cas_mismatch_cid",
    "cid_mismatch_cas",
    "cid_not_found",
}

# How many CIDs of a name/structure search are looked at, and the similarity
# (difflib ratio on normalized names) above which a synonym counts as close.
_MAX_SEARCH_CIDS = 3
_CLOSE_SYNONYM_RATIO = 0.8

# Stereodescriptors: "(R)", "(-)", "(+/-)", "(1S,5R)", "(E,E)", "cis-"...
_STEREO_GROUP_RE = re.compile(r"\(([0-9rsez,+\-±/ ]+)\)", re.IGNORECASE)
_CIS_TRANS_RE = re.compile(r"\b(cis|trans)\b", re.IGNORECASE)


def stereo_markers(name: str) -> set[str]:
    """The stereodescriptors a molecule name specifies, normalized so that
    equivalent spellings compare equal ("(±)" == "(+-)" == "(+/-)",
    "trans" == "E", "cis" == "Z")."""
    markers = {
        group.lower().replace(" ", "").replace("±", "+-").replace("+/-", "+-")
        for group in _STEREO_GROUP_RE.findall(name)
    }
    markers |= {"e" if m.lower() == "trans" else "z" for m in _CIS_TRANS_RE.findall(name)}
    return markers


@dataclass
class Row:
    """A Molecule row's identifiers, as typed in the Excel file."""

    name: str | None = None
    cids: list[int] = field(default_factory=list[int])
    cas: list[str] = field(default_factory=list[str])
    smiles: str | None = None
    inchikey: str | None = None
    mixture: str | None = None
    # The raw CAS cell: the pipeline resolves an empty-CID row's CID from it.
    cas_cell: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "cids": self.cids,
            "cas": self.cas,
            "smiles": self.smiles,
            "inchikey": self.inchikey,
            "mixture": self.mixture,
        }


@dataclass
class _Maps:
    names: dict[int, str]
    cas: dict[int, str]
    inchikeys: dict[int, str]
    smiles: dict[int, str]
    synonyms: dict[int, list[str]]

    @classmethod
    def for_cids(cls, cids: list[int]) -> _Maps:
        fetch_pubchem_data(sorted(set(cids)))
        return cls(*_parse_cid_cache(sorted(set(cids))))


def _finding(severity: str, code: str, message: str, **details: Any) -> dict[str, Any]:
    return {"severity": severity, "code": code, "message": message, **details}


def pipeline_outcome(row: Row) -> dict[str, Any]:
    """How main.py's molecule stage treats this row, step by step (see
    main._process_molecules and process_molecules.process_molecules)."""
    findings: list[dict[str, Any]] = []
    cids = list(row.cids)
    cas_list = list(row.cas)
    # validate_cid_or_cas counts any non-empty CAS cell, even unparseable.
    has_cas = bool(row.cas_cell) if row.cas_cell is not None else bool(cas_list)

    if not cids and not has_cas:
        if row.mixture == "mixture" and row.name:
            findings.append(
                _finding(
                    "warning",
                    "import_required_value_missing",
                    "Mixture row without CID/CAS: accepted by the pipeline, "
                    "skipped by the Laravel import (empty CID).",
                )
            )
        else:
            findings.append(
                _finding(
                    "error",
                    "cid_and_cas_missing",
                    "Molecule has neither a usable CID nor a CAS number.",
                )
            )
        return _summarize(findings, cids, row, None)

    if not cids:
        # main._process_molecules -> fetch_cids_from_cas: the whole CAS cell
        # is resolved to a CID.
        cas_key = normalize_dashes(row.cas_cell or ", ".join(cas_list))
        resolved = fetch_cas_to_cid_map([cas_key]).get(cas_key)
        if resolved:
            cids = [resolved]
            findings.append(
                _finding(
                    "info",
                    "cid_filled_from_cas",
                    "Empty CID filled from the CAS number.",
                    cid=resolved,
                )
            )
        else:
            findings.append(
                _finding(
                    "warning",
                    "import_required_value_missing",
                    "CAS doesn't resolve to a PubChem CID: the CID stays empty and "
                    "the Laravel import skips this row.",
                )
            )
            return _summarize(findings, cids, row, None)

    failed = set(fetch_pubchem_data(cids))
    if failed:
        findings.append(
            _finding(
                "warning",
                "cid_not_found",
                "CID not found on PubChem.",
                cids=sorted(failed),
            )
        )

    # CAS cross-check: single-CID rows only (get_unique_cas_for_cross_check).
    single = len(cids) == 1
    cas_to_cid = fetch_cas_to_cid_map(cas_list) if single else {}
    cas_details = fetch_cas_common_chemistry_details(cas_list) if single else {}
    maps = _Maps.for_cids([*cids, *cas_to_cid.values()])

    blocked = False
    cas_written = None
    if single and cas_list:
        check = check_cid_cas(
            row.name,
            cids[0],
            cas_list,
            cas_to_cid,
            maps.names,
            maps.synonyms,
            cas_details,
        )
        if check.outcome == "fix_cas":
            cas_written = maps.cas.get(cids[0])
            findings.append(
                _finding(
                    "auto_fix",
                    "cas_mismatch_cid",
                    "CAS resolves to a different compound; the Name matches the "
                    "CID, so the pipeline replaces the CAS with the CID's own.",
                    cas_to_cid=check.inconsistent,
                    pipeline_writes=cas_written,
                )
            )
        elif check.outcome == "fix_cid":
            new_cid = check.inconsistent[str(check.matching_cas)]
            findings.append(
                _finding(
                    "auto_fix",
                    "cid_mismatch_cas",
                    "CID doesn't match the CAS; the Name matches the CAS side, so "
                    "the pipeline replaces the CID.",
                    cas_to_cid=check.inconsistent,
                    pipeline_writes=new_cid,
                )
            )
            cids = [new_cid]
        elif check.outcome == "conflict":
            blocked = True
            findings.append(
                _finding(
                    "error",
                    "cid_cas_name_conflict",
                    "CID and CAS refer to different compounds and the Name matches "
                    "neither.",
                    cas_to_cid=check.inconsistent,
                )
            )

    if row.name and not blocked:
        matches = molecule_name_matches(
            row.name, cids, cas_list, cas_details, maps.names, maps.synonyms
        )
        if matches is False:
            findings.append(
                _finding(
                    "error",
                    "molecule_name_mismatch",
                    "Name not found among the PubChem name/synonyms of the CID(s), "
                    "nor the CAS Common Chemistry name/synonyms of the CAS number(s).",
                )
            )

    return _summarize(findings, cids, row, maps, cas_written)


def _summarize(
    findings: list[dict[str, Any]],
    cids: list[int],
    row: Row,
    maps: _Maps | None,
    cas_written: str | None = None,
) -> dict[str, Any]:
    status = "error" if any(f["severity"] == "error" for f in findings) else "passes"
    result: dict[str, Any] = {"status": status, "findings": findings, "final_cids": cids}
    if status == "passes" and maps is not None and cids:
        # What the export writes into the CSV (_fill_molecule_names,
        # _enrich_molecule_columns).
        written_cas = sorted({maps.cas[c] for c in cids if c in maps.cas})
        cas_details = (
            fetch_cas_common_chemistry_details(written_cas) if written_cas else {}
        )
        canonical = [
            canonical_name(c, maps.names, maps.cas, cas_details)[0] for c in cids
        ]
        spelled = (
            registry_spelling(
                row.name,
                cids[0],
                canonical[0],
                maps.names,
                maps.synonyms,
                maps.cas,
                cas_details,
            )
            if row.name and len(cids) == 1
            else None
        )
        result["pipeline_writes"] = {
            "Molecule Name": spelled or row.name,
            "Canonical Name": ", ".join(n for n in canonical if n)
            if all(canonical)
            else row.name,
            "CAS": cas_written
            or ", ".join(maps.cas[c] for c in cids if c in maps.cas)
            or None,
            "InChIKey": ", ".join(maps.inchikeys[c] for c in cids if c in maps.inchikeys)
            or None,
            "SMILES": ", ".join(maps.smiles[c] for c in cids if c in maps.smiles) or None,
        }
    return result


def _closest_synonyms(name: str, cid: int, maps: _Maps) -> list[dict[str, Any]]:
    target = _normalize_name(name)
    candidates = [maps.names[cid]] if cid in maps.names else []
    candidates += maps.synonyms.get(cid, [])
    scored = {
        synonym: difflib.SequenceMatcher(None, target, _normalize_name(synonym)).ratio()
        for synonym in candidates
    }
    best = sorted(scored.items(), key=lambda item: -item[1])[:3]
    return [{"synonym": s, "ratio": round(r, 3)} for s, r in best]


def _compound(cid: int, maps: _Maps, name: str | None) -> dict[str, Any]:
    smiles = maps.smiles.get(cid)
    stereo = analyze_stereo(smiles).status if smiles else None
    info: dict[str, Any] = {
        "cid": cid,
        "url": _PUBCHEM_COMPOUND_URL.format(cid=cid),
        "pubchem_title": maps.names.get(cid),
        "pubchem_cas": maps.cas.get(cid),
        "inchikey": maps.inchikeys.get(cid),
        "smiles": smiles,
        "stereo": stereo,
    }
    if name:
        info["name_matches"] = molecule_name_matches(
            name, [cid], [], {}, maps.names, maps.synonyms
        )
        info["closest_synonyms"] = _closest_synonyms(name, cid, maps)
    return info


def _lookups(row: Row) -> dict[str, list[int] | None]:
    lookups: dict[str, list[int] | None] = {}
    if row.name:
        lookups["name"] = fetch_cids_by_name(row.name)
    if row.smiles:
        lookups["smiles"] = fetch_cids_by_smiles(row.smiles)
    if row.inchikey:
        lookups["inchikey"] = fetch_cids_by_inchikey(row.inchikey)
    return lookups


def _suggestion(
    row: Row, changes: dict[str, Any], reason: str, evidence: list[str]
) -> dict[str, Any]:
    """A candidate correction of the row, re-checked with the pipeline rules."""
    fixed = Row(
        name=changes.get("Molecule Name", row.name),
        cids=[int(changes["CID"])] if "CID" in changes else row.cids,
        cas=[changes["CAS"]] if "CAS" in changes else row.cas,
        mixture=row.mixture,
    )
    before = {
        "Molecule Name": row.name,
        "CID": ", ".join(str(c) for c in row.cids) or None,
        "CAS": ", ".join(row.cas) or None,
    }
    return {
        "changes": [
            {"column": column, "from": before[column], "to": value}
            for column, value in changes.items()
            if before[column] != (str(value) if value is not None else None)
        ],
        "reason": reason,
        "evidence": evidence,
        "pipeline_after": pipeline_outcome(fixed)["status"],
    }


def _suggestions(
    row: Row, pipeline: dict[str, Any], lookups: dict[str, list[int] | None], maps: _Maps
) -> list[dict[str, Any]]:
    suggestions: list[dict[str, Any]] = []
    codes = {f["code"] for f in pipeline["findings"]}
    current = row.cids[0] if len(row.cids) == 1 else None

    for finding in pipeline["findings"]:
        if finding["code"] == "cas_mismatch_cid" and finding.get("pipeline_writes"):
            suggestions.append(
                _suggestion(
                    row,
                    {"CAS": finding["pipeline_writes"]},
                    "The CAS points to another compound; the Name confirms the CID, "
                    "whose PubChem record lists this CAS.",
                    [_PUBCHEM_COMPOUND_URL.format(cid=current)],
                )
            )
        if finding["code"] == "cid_mismatch_cas":
            new_cid = finding["pipeline_writes"]
            suggestions.append(
                _suggestion(
                    row,
                    {"CID": new_cid},
                    "The CID points to another compound; the Name confirms the CAS, "
                    "which resolves to this CID.",
                    [_PUBCHEM_COMPOUND_URL.format(cid=new_cid)],
                )
            )

    if not codes & {"molecule_name_mismatch", "cid_cas_name_conflict", "cid_and_cas_missing"}:
        return suggestions

    # The Name identifies another compound: take its CID (and PubChem's CAS
    # for it, when the row has a CAS column value to fix).
    found = (lookups.get("name") or [])[:_MAX_SEARCH_CIDS]
    for cid in found:
        if cid == current:
            continue
        changes: dict[str, Any] = {"CID": cid}
        if maps.cas.get(cid) and (row.cas or "cid_and_cas_missing" in codes):
            changes["CAS"] = maps.cas[cid]
        suggestions.append(
            _suggestion(
                row,
                changes,
                f"PubChem's name search resolves the Molecule Name to CID {cid} "
                f"({maps.names.get(cid)}).",
                [_PUBCHEM_COMPOUND_URL.format(cid=cid)]
                + (
                    [_CAS_COMMON_CHEMISTRY_DETAIL_URL.format(cas=changes["CAS"])]
                    if changes.get("CAS")
                    else []
                ),
            )
        )

    if current is not None and current in found:
        # PubChem itself says the name is this compound: the identifiers are
        # right and only the spelling trips the pipeline's name comparison.
        suggestions.append(
            _suggestion(
                row,
                {"Molecule Name": maps.names.get(current)},
                "PubChem's name search resolves the Molecule Name to this very CID, "
                "but the name isn't among its synonyms after normalization: the "
                "identifiers are consistent, only the spelling differs (possibly a "
                "gap in the pipeline's name normalization).",
                [_PUBCHEM_COMPOUND_URL.format(cid=current)],
            )
        )
    elif current is not None and row.name:
        # Probable typo: a synonym of the row's CID spelled almost the same -
        # never one naming another stereoisomer (e.g. "(R)-(-)-piperitone"
        # vs "(+-)-Piperitone"): that would change what was tested.
        for close in _closest_synonyms(row.name, current, maps):
            if close["ratio"] < _CLOSE_SYNONYM_RATIO:
                break
            if stereo_markers(close["synonym"]) != stereo_markers(row.name):
                continue
            suggestions.append(
                _suggestion(
                    row,
                    {"Molecule Name": close["synonym"]},
                    f"Close spelling (similarity {close['ratio']}) of a PubChem "
                    f"synonym of CID {current}: possible typo in the Name.",
                    [_PUBCHEM_COMPOUND_URL.format(cid=current)],
                )
            )
            break

    # Passing candidates first.
    return sorted(suggestions, key=lambda s: s["pipeline_after"] != "passes")


def lookup(row: Row) -> dict[str, Any]:
    collector = IssueCollector(catch_unexpected=False)
    with collector.stage("lookup"):
        pipeline = pipeline_outcome(row)
        lookups = _lookups(row)
        extra_cids = [c for cids in lookups.values() for c in (cids or [])[:_MAX_SEARCH_CIDS]]
        cas_to_cid = fetch_cas_to_cid_map(row.cas)
        cas_details = fetch_cas_common_chemistry_details(row.cas)
        all_cids = list(
            dict.fromkeys([*row.cids, *pipeline["final_cids"], *cas_to_cid.values(), *extra_cids])
        )
        maps = _Maps.for_cids(all_cids)
        suggestions = _suggestions(row, pipeline, lookups, maps)

    return _jsonable(
        {
            "input": row.to_json(),
            "pipeline": pipeline,
            "compounds": [_compound(cid, maps, row.name) for cid in all_cids],
            "cas": [
                {
                    "cas": cas,
                    "url": _CAS_COMMON_CHEMISTRY_DETAIL_URL.format(cas=cas),
                    "resolves_to_cid": cas_to_cid.get(cas),
                    "cas_common_chemistry": cas_details.get(cas),
                }
                for cas in row.cas
            ],
            "lookups": lookups,
            "suggestions": suggestions,
            "warnings": list(dict.fromkeys(issue.message for issue in collector.issues)),
        }
    )


# --- Excel rows ---------------------------------------------------------------


def _cell(value: object) -> str | None:
    if value is None or pd.isna(value) or not str(value).strip():  # pyright: ignore[reportUnknownArgumentType]
        return None
    return str(value).strip()


def rows_from_excel(excel_path: Path, excel_rows: list[int]) -> list[dict[str, Any]]:
    """The Molecule identifiers of the given Excel rows (data sheet)."""
    raw = get_raw_data_from_excel_file(excel_path)
    layout = ExcelLayout.from_raw(raw)
    columns = [MOLECULE_NAME, CID, CAS, INCHIKEY, SMILES, MIXTURE]
    found: list[dict[str, Any]] = []
    for excel_row in excel_rows:
        idx = excel_row - EXCEL_FIRST_DATA_ROW
        if idx not in raw.index:
            msg = f"Row {excel_row} is outside the data rows of {excel_path}"
            raise SystemExit(msg)
        mol = raw.loc[idx, MOLECULE]
        cas_cell = _cell(mol[CAS])
        mixture = _cell(mol[MIXTURE])
        row = Row(
            name=_cell(mol[MOLECULE_NAME]),
            cids=parse_cid_cell(mol[CID]),
            cas=parse_cas_cell(mol[CAS]),
            smiles=_cell(mol[SMILES]),
            inchikey=_cell(mol[INCHIKEY]),
            mixture=mixture.lower() if mixture else None,
            cas_cell=cas_cell,
        )
        found.append(
            {
                "excel": {
                    "file": str(excel_path),
                    "sheet": layout.sheet,
                    "row": excel_row,
                    "cells": {
                        column: f"{layout.excel_column(MOLECULE, column)}{excel_row}"
                        for column in columns
                        if layout.excel_column(MOLECULE, column)
                    },
                    "values": {column: _jsonable(mol[column]) for column in columns},
                },
                "row": row,
            }
        )
    return found


def lookup_report(report_path: Path) -> list[dict[str, Any]]:
    """Diagnose every molecule issue of a `main.py --check --report json`
    report, one lookup per distinct Excel row triple."""
    report = json.loads(report_path.read_text(encoding="utf-8"))
    issues = [i for i in report["issues"] if i["code"] in MOLECULE_ISSUE_CODES]
    first_rows = [i["excel_rows"][0] for i in issues if i["excel_rows"]]
    rows = {
        item["excel"]["row"]: item
        for item in rows_from_excel(Path(report["excel_path"]), first_rows)
    }
    results: list[dict[str, Any]] = []
    for issue in issues:
        if not issue["excel_rows"]:
            continue
        item = rows[issue["excel_rows"][0]]
        results.append(
            {
                "issue": {
                    key: issue[key]
                    for key in ("severity", "code", "column", "value", "cells", "excel_rows")
                },
                "excel": item["excel"],
                "lookup": lookup(item["row"]),
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--name", help="Molecule name")
    parser.add_argument("--cid", help="PubChem CID(s), as in a CID cell")
    parser.add_argument("--cas", help="CAS number(s), as in a CAS cell")
    parser.add_argument("--smiles", help="SMILES")
    parser.add_argument("--inchikey", help="InChIKey")
    parser.add_argument("--excel", type=Path, help="Study Excel file (with --row)")
    parser.add_argument(
        "--row", type=int, action="append", help="Excel row number (repeatable)"
    )
    parser.add_argument(
        "--report",
        type=Path,
        help="Check report (main.py --check --report json): diagnose all its "
        "molecule issues",
    )
    args = parser.parse_args()

    if args.report:
        result: Any = lookup_report(args.report)
    elif args.excel:
        if not args.row:
            parser.error("--excel requires at least one --row")
        result = [
            {"excel": item["excel"], "lookup": lookup(item["row"])}
            for item in rows_from_excel(args.excel, args.row)
        ]
    elif any((args.name, args.cid, args.cas, args.smiles, args.inchikey)):
        result = lookup(
            Row(
                name=args.name,
                cids=parse_cid_cell(args.cid) if args.cid else [],
                cas=parse_cas_cell(args.cas) if args.cas else [],
                smiles=args.smiles,
                inchikey=args.inchikey,
                cas_cell=args.cas,
            )
        )
    else:
        parser.error("give --name/--cid/--cas/--smiles/--inchikey, --excel/--row or --report")

    json.dump(result, sys.stdout, indent=2, ensure_ascii=False)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
