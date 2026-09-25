from typing import Any

import pytest

from scripts import cache_manager
from scripts.cache_manager import set_cache
from scripts.process_molecules import check_cid_cas, molecule_name_matches
from scripts.tools import lookup_molecule
from scripts.tools.lookup_molecule import Row, lookup, stereo_markers


def _record(title: str, cas: str, synonyms: list[str]) -> dict[str, Any]:
    """Minimal PubChem cache entry, as fetch_pubchem_data stores it."""
    return {
        "PC_Compounds": [
            {
                "props": [
                    {"urn": {"label": "InChIKey"}, "value": {"sval": f"KEY-{title}"}},
                    {
                        "urn": {"label": "SMILES", "name": "Absolute"},
                        "value": {"sval": "CC1CCC(CC1)O"},
                    },
                ]
            }
        ],
        "synonyms": [title, cas, *synonyms],
        "record_title": title,
    }


@pytest.fixture
def offline(tmp_path, monkeypatch):  # noqa: ANN001, ANN201
    """Fake PubChem cache and registries; no network call."""
    monkeypatch.setattr(cache_manager, "CACHE_DIR", str(tmp_path))
    set_cache("87468080", _record("Methylcyclohexanol", "25639-42-3", []), "molecules")
    set_cache(
        "11524",
        _record("4-Methylcyclohexanol", "589-91-3", ["p-methylcyclohexanol"]),
        "molecules",
    )
    cas_to_cid = {"25639-42-3": 87468080, "589-91-3": 11524}
    names = {"4-methylcyclohexanol": [11524]}
    monkeypatch.setattr(lookup_molecule, "fetch_pubchem_data", lambda cids: [])
    monkeypatch.setattr(
        lookup_molecule,
        "fetch_cas_to_cid_map",
        lambda cas: {c: cas_to_cid[c] for c in cas if c in cas_to_cid},
    )
    monkeypatch.setattr(
        lookup_molecule, "fetch_cas_common_chemistry_details", lambda cas: {}
    )
    monkeypatch.setattr(
        lookup_molecule, "fetch_cids_by_name", lambda name: names.get(name.lower(), [])
    )


def test_wrong_identifiers_get_a_passing_suggestion(offline):  # noqa: ANN001
    result = lookup(Row(name="4-Methylcyclohexanol", cas=["25639-42-3"]))
    codes = [f["code"] for f in result["pipeline"]["findings"]]
    assert result["pipeline"]["status"] == "error"
    assert codes == ["cid_filled_from_cas", "molecule_name_mismatch"]
    best = result["suggestions"][0]
    assert best["pipeline_after"] == "passes"
    assert {(c["column"], c["to"]) for c in best["changes"]} == {
        ("CID", 11524),
        ("CAS", "589-91-3"),
    }


def test_consistent_row_passes_without_suggestion(offline):  # noqa: ANN001
    result = lookup(Row(name="p-Methylcyclohexanol", cids=[11524], cas=["589-91-3"]))
    assert result["pipeline"]["status"] == "passes"
    assert result["pipeline"]["pipeline_writes"]["Molecule Name"] == "4-Methylcyclohexanol"
    assert result["suggestions"] == []


def test_typo_suggestion_keeps_identifiers(offline):  # noqa: ANN001
    result = lookup(Row(name="4-Methylcyclohexanl", cids=[11524]))
    (suggestion,) = result["suggestions"]
    assert suggestion["changes"] == [
        {"column": "Molecule Name", "from": "4-Methylcyclohexanl", "to": "4-Methylcyclohexanol"}
    ]
    assert suggestion["pipeline_after"] == "passes"


def test_cas_pointing_elsewhere_is_what_the_pipeline_rewrites(offline):  # noqa: ANN001
    result = lookup(Row(name="4-Methylcyclohexanol", cids=[11524], cas=["25639-42-3"]))
    (finding,) = result["pipeline"]["findings"]
    assert (finding["code"], finding["pipeline_writes"]) == ("cas_mismatch_cid", "589-91-3")
    assert result["suggestions"][0]["changes"] == [
        {"column": "CAS", "from": "25639-42-3", "to": "589-91-3"}
    ]


@pytest.mark.parametrize(
    ("name", "markers"),
    [
        ("(R)-(-)-piperitone", {"r", "-"}),
        ("(±)-Linalool", {"+-"}),
        ("(+/-)-Linalool", {"+-"}),
        ("trans-2-Hexenal", {"e"}),
        ("(E,E)-Farnesol", {"e,e"}),
        ("2-Ethyltoluene", set()),
    ],
)
def test_stereo_markers(name: str, markers: set[str]):
    assert stereo_markers(name) == markers


# --- the pipeline rules the tool shares -----------------------------------


_NAMES = {1: "Citronellal", 2: "(S)-(-)-Citronellal"}
_SYNONYMS = {1: ["citronellal", "106-23-0"], 2: ["(S)-citronellal", "5949-05-3"]}


@pytest.mark.parametrize(
    ("name", "cid", "outcome"),
    [
        ("Citronellal", 1, "fix_cas"),  # Name confirms the CID: the CAS is wrong
        ("(S)-(-)-Citronellal", 1, "fix_cid"),  # Name confirms the CAS side
        ("Geraniol", 1, "conflict"),  # Name confirms neither
        ("(S)-(-)-Citronellal", 2, "agree"),  # CAS resolves to the row's CID
    ],
)
def test_check_cid_cas(name: str, cid: int, outcome: str):
    check = check_cid_cas(name, cid, ["5949-05-3"], {"5949-05-3": 2}, _NAMES, _SYNONYMS, {})
    assert check.outcome == outcome


def test_molecule_name_matches_positionally_for_mixtures():
    names = {1: "Citronellal", 2: "Geraniol"}
    synonyms: dict[int, list[str]] = {}
    assert molecule_name_matches("Citronellal and Geraniol", [1, 2], [], {}, names, synonyms)
    assert not molecule_name_matches("Geraniol and Citronellal", [1, 2], [], {}, names, synonyms)
    assert molecule_name_matches("anything", [3], [], {}, names, synonyms) is None
