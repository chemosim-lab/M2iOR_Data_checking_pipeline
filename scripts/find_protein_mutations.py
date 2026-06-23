from Bio import Align

SUBSTITUTION_MATRIX = "BLOSUM62"


def _make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"  # Needleman-Wunsch
    aligner.substitution_matrix = Align.substitution_matrices.load(SUBSTITUTION_MATRIX)
    aligner.open_gap_score = -10
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


def align_and_annotate(seq: str, seq_ref: str) -> tuple[str, float, float]:
    """
    Align seq (query, from the study) against seq_ref (reference, from UniProt).

    Returns:
        mutations_str: mutations joined by "_", empty string if sequences are identical
        pid_aln: identity over alignment length (gaps included), rounded to 2 decimals
        pid_short: identity over the shorter sequence length, rounded to 2 decimals
    """
    aligner = _make_aligner()
    aln = aligner.align(seq_ref, seq)[0]
    a, b = aln[0], aln[1]

    mutations: list[str] = []
    ref_pos = 0
    for ref_aa, query_aa in zip(a, b, strict=True):
        if ref_aa != "-":
            ref_pos += 1
        if ref_aa == "-":
            mutations.append(f"ins{ref_pos}{query_aa}")
        elif query_aa == "-":
            mutations.append(f"{ref_aa}{ref_pos}del")
        elif ref_aa != query_aa:
            mutations.append(f"{ref_aa}{ref_pos}{query_aa}")

    matches = sum(x == y and x != "-" for x, y in zip(a, b, strict=True))
    aln_len = len(a)
    shorter = min(len(seq), len(seq_ref))
    pid_aln = round(100 * matches / aln_len, 2) if aln_len > 0 else 0.0
    pid_short = round(100 * matches / shorter, 2) if shorter > 0 else 0.0

    return "_".join(mutations), pid_aln, pid_short


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
