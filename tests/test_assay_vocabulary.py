"""Allowed values of the free-text Assay columns (scripts/process_assay.py)."""

import pandas as pd
import pytest

from scripts import process_assay
from scripts.columns import (
    ASSAY,
    EXPERIMENTAL_TECHNIQUE,
    EXPRESSION_DRIVER,
    EXPRESSION_SYSTEM,
    RECORDING_SYSTEM,
    TYPE,
)
from scripts.normalize_dataframe import normalize_df
from scripts.process_assay import (
    load_vocabulary,
    normalize_assay_field_format,
    parse_vocabulary,
    registry_key,
    validate_assay_consistency,
    validate_assay_vocabulary,
)
from scripts.process_responses import (
    _ALLOWED_EXPERIMENTAL_TECHNIQUES,
    validate_main_flux_unit_column,
    validate_stimulation_duration_unit_column,
)
from scripts.report import IssueCollector, ValidationError

SSR = "single sensillum recording"
TEVC = "two-electrode voltage clamp"

_TEST_VOCABULARY = {
    RECORDING_SYSTEM: {
        "definition": "Device that records the response.",
        "allowed": {
            "OC-725C": {"requires": {EXPERIMENTAL_TECHNIQUE: [TEVC]}},
            "tungsten electrode": {"requires": {EXPERIMENTAL_TECHNIQUE: [SSR]}},
        },
        "pending": {"Manual": "Not in any article."},
    },
    TYPE: {
        "allowed": {"∆ fluorescence": {}, "spiking frequency": {}},
    },
    EXPRESSION_DRIVER: {
        "allowed": {
            "Or22a-Gal4": {"requires": {EXPRESSION_SYSTEM: ["Drosophila"]}},
            "pT7Ts": {"requires": {EXPRESSION_SYSTEM: ["Xenopus oocytes"]}},
        },
    },
}


@pytest.fixture(autouse=True)
def _vocabulary(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        process_assay, "_VOCABULARY", parse_vocabulary(_TEST_VOCABULARY)
    )


def _frame(columns: dict[str, list[object]]) -> pd.DataFrame:
    return pd.DataFrame({(ASSAY, name): values for name, values in columns.items()})


def _collect(step, *args) -> IssueCollector:  # noqa: ANN001
    collector = IssueCollector(catch_unexpected=True)
    with collector.stage("test"):
        step(*args)
    return collector


def test_registry_key_ignores_case_whitespace_and_delta_sign():
    assert registry_key("  OC-725C  Amplifier ") == "oc-725c amplifier"
    assert registry_key("Δ fluorescence") == registry_key("∆ fluorescence")


def test_case_only_difference_is_respelled_once_per_value():
    df = _frame({RECORDING_SYSTEM: ["oc-725c", "OC-725c", "OC-725C", None]})
    collector = _collect(normalize_assay_field_format, df)

    assert df[ASSAY][RECORDING_SYSTEM].tolist()[:3] == ["OC-725C"] * 3
    assert pd.isna(df[ASSAY][RECORDING_SYSTEM].iloc[3])
    issues = {i.value: i for i in collector.issues}
    assert set(issues) == {"oc-725c", "OC-725c"}
    assert issues["oc-725c"].code == "assay_field_format"
    assert issues["oc-725c"].severity == "auto_fix"
    assert issues["oc-725c"].suggested_value == "OC-725C"
    assert issues["oc-725c"].rows == [0]


def test_greek_delta_is_respelled_to_the_listed_increment_sign():
    df = _frame({TYPE: ["δ fluorescence"]})
    _collect(normalize_assay_field_format, df)
    assert df[ASSAY][TYPE].tolist() == ["∆ fluorescence"]


def test_pending_value_is_kept_and_only_warns():
    df = _frame({RECORDING_SYSTEM: ["Manual"]})
    collector = _collect(normalize_assay_field_format, df)

    assert df[ASSAY][RECORDING_SYSTEM].tolist() == ["Manual"]
    (issue,) = collector.issues
    assert issue.code == "assay_value_pending"
    assert issue.severity == "warning"
    validate_assay_vocabulary(df)  # does not block


def test_value_outside_the_list_blocks_with_the_closest_listed_value():
    df = _frame({RECORDING_SYSTEM: ["OC-725C", "OC-725 C amplifier", "OC-7256"]})
    with pytest.raises(ValidationError) as caught:
        validate_assay_vocabulary(df)

    issues = {i.value: i for i in caught.value.issues}
    assert set(issues) == {"OC-725 C amplifier", "OC-7256"}
    issue = issues["OC-7256"]
    assert issue.code == "value_not_allowed"
    assert issue.severity == "error"
    assert issue.suggested_value == "OC-725C"
    assert issue.details["allowed"] == ["OC-725C", "tungsten electrode"]
    assert issue.rows == [2]


def test_empty_cells_are_not_checked():
    df = _frame({RECORDING_SYSTEM: [None, "", "  "]})
    validate_assay_vocabulary(df)


def test_value_on_a_row_of_another_technique_is_inconsistent():
    df = _frame(
        {
            EXPERIMENTAL_TECHNIQUE: [TEVC, "calcium imaging", SSR],
            RECORDING_SYSTEM: ["OC-725C", "OC-725C", "tungsten electrode"],
        }
    )
    with pytest.raises(ValidationError) as caught:
        validate_assay_consistency(df)

    (issue,) = caught.value.issues
    assert issue.code == "assay_inconsistent"
    assert issue.column == RECORDING_SYSTEM
    assert issue.value == "OC-725C"
    assert issue.rows == [1]
    assert issue.details["found"] == "calcium imaging"


def test_driver_must_match_the_expression_system():
    df = _frame(
        {
            EXPRESSION_SYSTEM: ["Drosophila", "Drosophila", "Xenopus oocytes"],
            EXPRESSION_DRIVER: ["Or22a-Gal4", "pT7Ts", "pT7Ts"],
        }
    )
    with pytest.raises(ValidationError) as caught:
        validate_assay_consistency(df)

    (issue,) = caught.value.issues
    assert issue.value == "pT7Ts"
    assert issue.rows == [1]


def test_pending_values_are_not_checked_for_consistency():
    df = _frame(
        {
            EXPERIMENTAL_TECHNIQUE: ["calcium imaging"],
            RECORDING_SYSTEM: ["Manual"],
        }
    )
    validate_assay_consistency(df)


# --- the shipped vocabulary ---------------------------------------------------


def test_shipped_vocabulary_is_self_consistent():
    vocabulary = load_vocabulary()
    systems = {
        e.canonical for e in vocabulary[EXPRESSION_SYSTEM].entries.values()
    }
    for column, column_vocabulary in vocabulary.items():
        assert column_vocabulary.definition, column
        keys = [registry_key(e.canonical) for e in column_vocabulary.entries.values()]
        assert len(keys) == len(set(keys)), f"{column}: values differing only in case"
        for entry in column_vocabulary.entries.values():
            if entry.pending:
                assert entry.note, f"{column}: pending {entry.canonical!r} needs a note"
                assert not entry.requires
            for other, expected in entry.requires.items():
                assert other in (EXPERIMENTAL_TECHNIQUE, EXPRESSION_SYSTEM)
                if other == EXPERIMENTAL_TECHNIQUE:
                    assert set(expected) <= _ALLOWED_EXPERIMENTAL_TECHNIQUES
                else:
                    assert set(expected) <= systems


def test_shipped_vocabulary_is_stored_in_its_own_normalized_spelling():
    """Lowercased columns (Type, Odor delivery system) must be listed in lowercase:
    cells are lowercased before they are matched."""
    vocabulary = load_vocabulary()
    for column in (TYPE, "Odor delivery system"):
        for entry in vocabulary[column].entries.values():
            assert entry.canonical == entry.canonical.lower()


# --- unit lists ---------------------------------------------------------------


def test_second_is_an_alias_of_s_in_the_duration_unit():
    df = pd.DataFrame({(ASSAY, "Stimulation duration unit"): ["second", "ms", "s"]})
    df.columns = pd.MultiIndex.from_tuples(df.columns)
    df = normalize_df(df)
    assert df[ASSAY]["Stimulation duration unit"].tolist() == ["s", "ms", "s"]
    validate_stimulation_duration_unit_column(df)


def test_main_flux_unit_is_validated_like_the_stimulation_flux_unit():
    df = pd.DataFrame({(ASSAY, "Main Flux Unit"): ["l/min", "gallons"]})
    df.columns = pd.MultiIndex.from_tuples(df.columns)
    with pytest.raises(ValidationError) as caught:
        validate_main_flux_unit_column(df)
    assert caught.value.issues[0].value == "gallons"
