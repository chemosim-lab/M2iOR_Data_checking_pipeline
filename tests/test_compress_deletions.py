import pytest

from scripts.find_protein_mutations import compress_deletions


def test_all_consecutive_deletions():
    assert compress_deletions(
        "F300del;Q301del;V302del;F303del;F304del;Y305del"
    ) == "F300_Y305del"


def test_substitutions_then_consecutive_deletions():
    assert compress_deletions(
        "R200H;T201A;G296del;C297del;M298del;L299del;F300del;Q301del;V302del;F303del"
    ) == "R200H;T201A;G296_F303del"


def test_deletions_with_position_gap():
    assert compress_deletions("G296del;Q301del") == "G296del;Q301del"


def test_gap_then_consecutive_deletions():
    assert compress_deletions("G296del;M298del;L299del") == "G296del;M298_L299del"


def test_substitution_breaks_deletion_run():
    assert compress_deletions("G296del;X297Y;M298del") == "G296del;X297Y;M298del"


def test_unknown_token_raises():
    with pytest.raises(ValueError, match="P100fs"):
        compress_deletions("P100fs")


def test_single_deletion_unchanged():
    assert compress_deletions("A42del") == "A42del"


def test_single_substitution_unchanged():
    assert compress_deletions("R200H") == "R200H"


def test_empty_string():
    assert compress_deletions("") == ""


def test_unknown_token_in_middle_raises():
    with pytest.raises(ValueError, match="K50dup"):
        compress_deletions("A42del;K50dup;L51del")


def test_two_separate_runs():
    assert compress_deletions("A1del;B2del;D5del;E6del;F7del") == "A1_B2del;D5_F7del"
