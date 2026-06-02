from Bio import Align


def _make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"  # Needleman-Wunsch
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    return aligner


def find_mutations(seq1: str, seq2: str, seq1_name="Ref", seq2_name="Query") -> list[str]:
    """
    Compare deux séquences protéiques et retourne la liste des mutations.
    Format de sortie : AncienAAPositionNouvelAA (ex: A42G)
    seq1 = reference, seq2 = query
    """
    aligner = _make_aligner()
    alignments = aligner.align(seq1, seq2)
    best = alignments[0]

    mutations: list[str] = []
    ref_pos = 0

    for ref_aa, query_aa in zip(*best):
        if ref_aa != "-":
            ref_pos += 1
        if ref_aa == "-":
            mutations.append(f"ins{ref_pos}{query_aa}")
        elif query_aa == "-":
            mutations.append(f"{ref_aa}{ref_pos}del")
        elif ref_aa != query_aa:
            mutations.append(f"{ref_aa}{ref_pos}{query_aa}")

    return mutations


def align_and_annotate(seq: str, seq_ref: str) -> tuple[str, float]:
    """
    Align seq (query, from the study) against seq_ref (reference, from UniProt).

    Returns:
        mutations_str: mutations joined by "_", empty string if sequences are identical
        identity_pct: identical positions / alignment length * 100, rounded to 2 decimals
    """
    aligner = _make_aligner()
    alignments = aligner.align(seq_ref, seq)
    best = alignments[0]

    mutations: list[str] = []
    ref_pos = 0
    identical = 0
    total = 0

    for ref_aa, query_aa in zip(*best):
        total += 1
        if ref_aa != "-":
            ref_pos += 1
        if ref_aa == "-":
            mutations.append(f"ins{ref_pos}{query_aa}")
        elif query_aa == "-":
            mutations.append(f"{ref_aa}{ref_pos}del")
        elif ref_aa != query_aa:
            mutations.append(f"{ref_aa}{ref_pos}{query_aa}")
        else:
            identical += 1

    identity = round(identical / total * 100, 2) if total > 0 else 0.0
    return "_".join(mutations), identity


if __name__ == "__main__":
    _ref = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGEDEDTLLCIDVDPQMEIHDGVPFTLEDFKPHHVKKFDELDLH"
    _query = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGEDEDTLLCIDVDPQMEIHDGVPFTLEDFKPHHVKKFDELDLH"
    _mutations = find_mutations(_ref, _query)
    if _mutations:
        print(f"{len(_mutations)} mutation(s) trouvée(s) :")
        for m in _mutations:
            print(f"  - {m}")
    else:
        print("Aucune mutation — séquences identiques.")
