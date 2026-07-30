import re

from Bio import Align

SUBSTITUTION_MATRIX = "BLOSUM62"

_DEL_RE = re.compile(r"^([A-Z])(\d+)del$")
_SUB_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
# HGVS: X10_Y11insABC  — fallback for terminal insertions: ins0ABC
_INS_RE = re.compile(r"^(?:[A-Z]\d+_[A-Z]\d+ins[A-Z]+|ins\d+[A-Z]+)$")


def compress_deletions(mutations_str: str) -> str:
    """
    Compress runs of consecutive single-residue deletions into HGVS range notation.

    "F300del;Q301del;V302del" -> "F300_V302del"

    Token separator is `;`. The `_` character is reserved for deletion ranges.
    Substitutions and insertions pass through unchanged and break any ongoing run.
    Any token that is neither a deletion, substitution, nor insertion raises ValueError.
    """
    if not mutations_str:
        return mutations_str

    tokens = mutations_str.split(";")
    output: list[str] = []
    run: list[tuple[str, int]] = []

    def _flush() -> None:
        if not run:
            return
        if len(run) == 1:
            aa, pos = run[0]
            output.append(f"{aa}{pos}del")
        else:
            aa0, pos0 = run[0]
            aa1, pos1 = run[-1]
            output.append(f"{aa0}{pos0}_{aa1}{pos1}del")
        run.clear()

    for token in tokens:
        del_m = _DEL_RE.match(token)
        if del_m:
            aa, pos = del_m.group(1), int(del_m.group(2))
            if run and pos != run[-1][1] + 1:
                _flush()
            run.append((aa, pos))
            continue

        if _SUB_RE.match(token) or _INS_RE.match(token):
            _flush()
            output.append(token)
            continue

        msg = f"Unrecognized mutation token: {token!r}"
        raise ValueError(msg)

    _flush()
    return ";".join(output)


def _emit_insertion(seq_ref: str, ref_pos: int, ins_aas: list[str]) -> str:
    """Build a HGVS protein insertion token for consecutive inserted residues."""
    inserted = "".join(ins_aas)
    # Terminal insertions have no flanking pair — use simplified fallback notation
    if ref_pos == 0 or ref_pos >= len(seq_ref):
        return f"ins{ref_pos}{inserted}"
    aa_left = seq_ref[ref_pos - 1]
    aa_right = seq_ref[ref_pos]
    return f"{aa_left}{ref_pos}_{aa_right}{ref_pos + 1}ins{inserted}"


def _make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Smith-Waterman
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
    ins_buffer: list[str] = []

    for ref_aa, query_aa in zip(*best):
        if ref_aa != "-":
            if ins_buffer:
                mutations.append(_emit_insertion(seq1, ref_pos, ins_buffer))
                ins_buffer.clear()
            ref_pos += 1
            if query_aa == "-":
                mutations.append(f"{ref_aa}{ref_pos}del")
            elif ref_aa != query_aa:
                mutations.append(f"{ref_aa}{ref_pos}{query_aa}")
        else:
            ins_buffer.append(query_aa)

    if ins_buffer:
        mutations.append(_emit_insertion(seq1, ref_pos, ins_buffer))

    return mutations


def align_and_annotate(seq: str, seq_ref: str) -> tuple[str, float, float]:
    """
    Align seq (query, from the study) against seq_ref (reference, from UniProt).

    Returns:
        mutations_str: mutations joined by ";", empty string if sequences are identical
        pid_aln: identity over alignment length (gaps included), rounded to 2 decimals
        pid_short: identity over the shorter sequence length, rounded to 2 decimals

    Residues of either sequence that fall outside the local alignment's window
    (terminal indels) are deliberately NOT reported as mutations: unlike a
    substitution or an internal indel, a truncated/extended terminus is as
    likely to come from a signal peptide, an expression tag, an isoform, or a
    partial construct as from an actual mutation, so folding it into the same
    HGVS-token field as real substitutions would be misleading. Its effect is
    still visible as a gap between pid_aln (identity within the aligned core)
    and pid_short (identity over the shorter full-length sequence).
    """
    # A trailing "*" is a stop-codon marker from translation, not a residue.
    seq = seq.rstrip("*")
    seq_ref = seq_ref.rstrip("*")

    aligner = _make_aligner()
    aln = aligner.align(seq_ref, seq)[0]
    a, b = aln[0], aln[1]

    # A local alignment doesn't have to start at seq_ref's first residue:
    # ref_pos must be seeded at the alignment's actual start, not 0, or every
    # position reported from here on is silently offset.
    ref_blocks, _query_blocks = aln.aligned
    ref_pos = int(ref_blocks[0][0])

    mutations: list[str] = []
    ins_buffer: list[str] = []

    for ref_aa, query_aa in zip(a, b, strict=True):
        if ref_aa != "-":
            if ins_buffer:
                mutations.append(_emit_insertion(seq_ref, ref_pos, ins_buffer))
                ins_buffer.clear()
            ref_pos += 1
            if query_aa == "-":
                mutations.append(f"{ref_aa}{ref_pos}del")
            elif ref_aa != query_aa:
                mutations.append(f"{ref_aa}{ref_pos}{query_aa}")
        else:
            ins_buffer.append(query_aa)

    if ins_buffer:
        mutations.append(_emit_insertion(seq_ref, ref_pos, ins_buffer))

    matches = sum(x == y and x != "-" for x, y in zip(a, b, strict=True))
    aln_len = len(a)
    shorter = min(len(seq), len(seq_ref))
    pid_aln = round(100 * matches / aln_len, 2) if aln_len > 0 else 0.0
    pid_short = round(100 * matches / shorter, 2) if shorter > 0 else 0.0

    return compress_deletions(";".join(mutations)), pid_aln, pid_short


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
