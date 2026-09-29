# pipeline/scripts/process_responses.py
import difflib
import math
import numbers
import re

import pandas as pd

from scripts.columns import (
    ASSAY,
    CONCENTRATION,
    CONCENTRATION_UNIT,
    EXPERIMENTAL_TECHNIQUE,
    PARAMETER,
    RESPONSE,
    RESPONSIVE,
    STIMULATION_DURATION_UNIT,
    STIMULATION_FLUX_UNIT,
    UNIT,
    VALUE,
    VALUE_NATURE,
)
from scripts.report import Issue, ValidationError

# Forms of a Value/Concentration cell the website import reads (resolveMinMaxFloat).
_UNICODE_DASHES = str.maketrans(dict.fromkeys("\u2010\u2011\u2012\u2013\u2014\u2212", "-"))
_PHP_NUMERIC = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
_RANGE = re.compile(r"-?\d+(?:\.\d+)?--?\d+(?:\.\d+)?")
_ACCEPTED_NUMBER_FORMS = "a number, '<X', '<=X', '>X', '>=X', '=X' or a range 'X1-X2'"

_ALLOWED_PARAMETERS = {"primary", "secondary", "dose-response"}
_ALLOWED_VALUE_NATURES = {"raw", "norm_rec", "norm_pair", "norm_other", "ec50"}
_ALLOWED_UNITS = {
    "na",
    "mol/l",
    "spikes/s",
    "v/v",
    "ug/ul",
    "um",
    "m",
    "percent",
    "%δ fluorescence",
    "nm",
}
_ALLOWED_CONCENTRATION_UNITS = {
    "v/v",
    "m",
    "mg/ml",
    "mol/l",
    "ng/ul",
    "um",
    "ug/ul",
    "ug",
    "mg",
    "nmol",
    "nm",
}
_ALLOWED_STIMULATION_FLUX_UNITS = {"ml/s", "l/min", "ml/min"}
_ALLOWED_STIMULATION_DURATION_UNITS = {"ms", "s", "second"}
_ALLOWED_EXPERIMENTAL_TECHNIQUES = {
    "two-electrode voltage clamp",
    "calcium imaging",
    "single sensillum recording",
    "fluorescence",
    "door 2.0",
    "electroantennography",
    "hek293",
}


def _validate_allowed_values(
    df: pd.DataFrame, group: str, col_name: str, allowed: set[str]
) -> None:
    col = df[group][col_name].dropna().astype(str).str.strip().str.lower()
    invalid = col[~col.isin(allowed)]
    if not invalid.empty:
        bad = invalid.unique().tolist()
        msg = (
            f"Column '{col_name}' contains invalid value(s) {bad} "
            f"(allowed: {sorted(allowed)})"
        )
        issues: list[Issue] = []
        for value in bad:
            closest = difflib.get_close_matches(value, sorted(allowed), n=1, cutoff=0.6)
            issues.append(
                Issue(
                    code="value_not_allowed",
                    message=f"'{value}' is not an allowed '{col_name}' value.",
                    group=group,
                    column=col_name,
                    rows=invalid.index[invalid == value].tolist(),
                    value=value,
                    suggested_value=closest[0] if closest else None,
                    details={"allowed": sorted(allowed)},
                )
            )
        raise ValidationError(msg, issues)


def _row_issues(
    invalid: pd.Series, code: str, message: str, group: str, col_name: str
) -> list[Issue]:
    return [
        Issue(
            code=code,
            message=message,
            group=group,
            column=col_name,
            rows=[idx],
            value=value,
        )
        for idx, value in invalid.items()
    ]


def validate_responsive_column(df: pd.DataFrame) -> None:
    """Raise ValueError if any non-null value in Responsive is not 0 or 1."""
    col = df[RESPONSE][RESPONSIVE].dropna()
    numeric: pd.Series = pd.to_numeric(col, errors="coerce")
    invalid = col[~numeric.isin([0, 1]) | numeric.isna()]
    if not invalid.empty:
        rows = invalid.index.tolist()
        msg = (
            f"Column '{RESPONSIVE}' contains invalid value(s) (expected 0 or 1) "
            f"at row(s): {rows}"
        )
        issues = _row_issues(
            invalid,
            "invalid_responsive",
            f"'{RESPONSIVE}' must be 0 or 1.",
            RESPONSE,
            RESPONSIVE,
        )
        raise ValidationError(msg, issues)


def _validate_float_column(df: pd.DataFrame, col_name: str) -> None:
    col = df[RESPONSE][col_name].dropna()
    numeric: pd.Series = pd.to_numeric(col, errors="coerce")
    invalid = col[numeric.isna()]
    if not invalid.empty:
        rows = invalid.index.tolist()
        raise ValueError(
            f"Column '{col_name}' contains non-float value(s) at row(s): {rows}"
        )


def import_readable_number(value: object) -> bool:
    """Whether the website's import can read `value` as a Value/Concentration cell.

    Mirrors CsvImportService::resolveMinMaxFloat() of the website: after removing
    whitespace and turning Unicode dashes into '-', a cell is empty, a number
    (PHP is_numeric), '<X', '<=X', '>X', '>=X', '=X' or a range 'X1-X2'. The
    import skips a row whose cell is anything else ('~12', 'N.A.', '1,5', a date).
    """
    if isinstance(value, bool):
        return False
    if isinstance(value, numbers.Real):
        return math.isfinite(value)
    text = re.sub(r"\s+", "", str(value).translate(_UNICODE_DASHES))
    if text[:1] in ("<", ">"):
        bound = text[2:] if text[1:2] == "=" else text[1:]
        return bound == "" or _PHP_NUMERIC.fullmatch(bound) is not None
    if text[:1] == "=":
        text = text[1:]
    return (
        text == ""
        or _RANGE.fullmatch(text) is not None
        or _PHP_NUMERIC.fullmatch(text) is not None
    )


def _validate_import_readable(df: pd.DataFrame, col_name: str) -> None:
    col = df[RESPONSE][col_name].dropna()
    invalid = col[~col.map(import_readable_number).astype(bool)]
    if invalid.empty:
        return
    rows = invalid.index.tolist()
    msg = f"Column '{col_name}' contains value(s) the import can't read at row(s): {rows}"
    issues = [
        Issue(
            code="invalid_value",
            message=f"The website import can't read this '{col_name}' and would skip "
            "the row: write a number, '<X', '>X' or a range 'X1-X2'.",
            group=RESPONSE,
            column=col_name,
            rows=invalid.index[invalid.astype(str) == str(value)].tolist(),
            value=value,
            details={"accepted": _ACCEPTED_NUMBER_FORMS},
        )
        for value in dict.fromkeys(invalid.tolist())
    ]
    raise ValidationError(msg, issues)


def validate_value_column(df: pd.DataFrame) -> None:
    """Raise ValidationError if a Value can't be read by the website import."""
    _validate_import_readable(df, VALUE)


def validate_ec50_values(df: pd.DataFrame) -> None:
    """Raise ValidationError if an EC50 value isn't positive.

    An EC50 is a concentration: a value <= 0 is usually its log10 (log EC50)
    or a placeholder, which the website would show as a concentration.
    """
    resp = df[RESPONSE]
    nature = resp[VALUE_NATURE].astype(str).str.strip().str.lower()
    value: pd.Series = pd.to_numeric(resp[VALUE], errors="coerce")
    invalid = resp[VALUE][(nature == "ec50") & (value <= 0)]
    if not invalid.empty:
        rows = invalid.index.tolist()
        msg = f"Column '{VALUE}' contains EC50 value(s) <= 0 at row(s): {rows}"
        issues = [
            Issue(
                code="ec50_not_positive",
                message="An EC50 is a concentration and must be > 0; a negative "
                "value is usually a log10 (log EC50), which needs converting.",
                group=RESPONSE,
                column=VALUE,
                rows=[idx],
                value=raw,
                details={"unit": None if pd.isna(unit) else unit},
            )
            for idx, raw, unit in zip(
                invalid.index, invalid, resp[UNIT][invalid.index], strict=True
            )
        ]
        raise ValidationError(msg, issues)


def validate_concentration_column(df: pd.DataFrame) -> None:
    """Raise ValidationError if a Concentration can't be read by the website import."""
    _validate_import_readable(df, CONCENTRATION)


def validate_parameter_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(df, RESPONSE, PARAMETER, _ALLOWED_PARAMETERS)


def validate_value_nature_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(df, RESPONSE, VALUE_NATURE, _ALLOWED_VALUE_NATURES)


def validate_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(df, RESPONSE, UNIT, _ALLOWED_UNITS)


def validate_concentration_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, RESPONSE, CONCENTRATION_UNIT, _ALLOWED_CONCENTRATION_UNITS
    )


def validate_stimulation_flux_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, ASSAY, STIMULATION_FLUX_UNIT, _ALLOWED_STIMULATION_FLUX_UNITS
    )


def validate_stimulation_duration_unit_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, ASSAY, STIMULATION_DURATION_UNIT, _ALLOWED_STIMULATION_DURATION_UNITS
    )


def validate_experimental_technique_column(df: pd.DataFrame) -> None:
    _validate_allowed_values(
        df, ASSAY, EXPERIMENTAL_TECHNIQUE, _ALLOWED_EXPERIMENTAL_TECHNIQUES
    )
