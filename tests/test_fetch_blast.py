import pytest

from scripts import fetch_blast
from scripts.fetch_blast import accession_from_seq_id, fetch_blast_reference


@pytest.mark.parametrize(
    ("seq_id", "accession"),
    [
        ("ref|XP_081832551.1|", "XP_081832551.1"),
        ("gb|QOI12088.1|", "QOI12088.1"),
        ("emb|CAB3514346.1|", "CAB3514346.1"),
        ("sp|P81909.1|OR22A_DROME", "P81909.1"),
        ("gi|1234567|ref|XP_019768125.2|", "XP_019768125.2"),
        ("XP_081832551.1", "XP_081832551.1"),
        (" A0A7L8XYY7 ", "A0A7L8XYY7"),
    ],
)
def test_accession_from_seq_id(seq_id, accession):
    assert accession_from_seq_id(seq_id) == accession


def test_cached_raw_identifier_is_returned_as_accession(monkeypatch):
    # Cache entries written before the fix hold the raw BLAST identifier.
    key = fetch_blast._blast_key("MKTAYIAK", "Ips typographus")
    cache = {key: {"accession": "ref|XP_081832551.1|", "sequence_ref": "MKTAYIAK"}}
    monkeypatch.setattr(fetch_blast, "_load_blast_cache", lambda: cache)

    assert fetch_blast_reference("MKTAYIAK", "Ips typographus") == (
        "XP_081832551.1",
        "MKTAYIAK",
    )


def test_cached_failure_returns_none(monkeypatch):
    key = fetch_blast._blast_key("MKTAYIAK", "Ips typographus")
    cache = {key: {"accession": None, "sequence_ref": None, "error": "no hit found"}}
    monkeypatch.setattr(fetch_blast, "_load_blast_cache", lambda: cache)

    assert fetch_blast_reference("MKTAYIAK", "Ips typographus") is None
