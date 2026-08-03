import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

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


def _blastp_path() -> str:
    """
    Locate the blastp binary: first the current venv's bin/ (where
    `scripts/install_blast.sh` puts it), falling back to PATH.
    """
    venv_blastp = Path(sys.prefix) / "bin" / "blastp"
    if venv_blastp.exists():
        return str(venv_blastp)
    found = shutil.which("blastp")
    if found:
        return found
    msg = (
        "blastp introuvable (ni dans .venv/bin, ni dans le PATH). "
        "Lancer scripts/install_blast.sh pour l'installer."
    )
    raise RuntimeError(msg)


# Biopython's PairwiseAligner computes the exact, globally-optimal local
# alignment (Gotoh DP, no early termination) — unlike BLAST, which stops
# extending once its X-drop heuristic sees the score fall too far below the
# running max, and applies composition-based matrix adjustment. That gap
# matters for chimeric/misassembled reference records: if two genuinely
# homologous blocks are separated by a stretch of unrelated sequence (e.g. a
# fused neighboring gene, or a missing exon), Biopython will happily bridge
# the whole thing with one huge gap because the extra matched residues on the
# far side outscore the gap penalty — producing one alignment instead of
# BLAST's several separate HSPs, and tanking %identity in the process.
# Reproducing BLAST's own X-drop/composition-adjusted extension by hand isn't
# practical, so shell out to the real blastp instead: exact same numbers as
# blast.ncbi.nlm.nih.gov, by construction.
def _best_hsp(seq_ref: str, seq: str) -> dict[str, str | int] | None:
    """
    Run blastp with seq as query and seq_ref as subject, and return the
    best-scoring HSP's fields (blastp orders HSPs by decreasing bit score,
    so the first output row is always the best one), or None if blastp finds
    no hit at all.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        query_fasta = Path(tmpdir) / "query.fasta"
        subject_fasta = Path(tmpdir) / "subject.fasta"
        query_fasta.write_text(f">query\n{seq}\n")
        subject_fasta.write_text(f">subject\n{seq_ref}\n")
        result = subprocess.run(
            [
                _blastp_path(),
                "-query",
                str(query_fasta),
                "-subject",
                str(subject_fasta),
                "-outfmt",
                "6 nident length sstart sseq qseq",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

    lines = result.stdout.strip().splitlines()
    if not lines:
        return None
    nident, length, sstart, sseq, qseq = lines[0].split("\t")
    return {
        "nident": int(nident),
        "length": int(length),
        "sstart": int(sstart),
        "sseq": sseq,
        "qseq": qseq,
    }


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
    Align seq (query, from the study) against seq_ref (reference, from UniProt),
    via blastp, so pid_aln matches what blast.ncbi.nlm.nih.gov itself reports.

    Returns:
        mutations_str: mutations joined by ";", empty string if sequences are identical
        pid_aln: identity over the best HSP's alignment length (gaps included),
            rounded to 2 decimals — this is BLAST's own "Per. Ident" for the hit
        pid_short: identity over the shorter sequence length, rounded to 2 decimals

    Residues of either sequence that fall outside the best HSP's window
    (terminal indels, or an unrelated region entirely — e.g. a chimeric or
    misassembled reference record) are deliberately NOT reported as
    mutations: unlike a substitution or an internal indel, a truncated/
    extended terminus is as likely to come from a signal peptide, an
    expression tag, an isoform, or a partial construct as from an actual
    mutation, so folding it into the same HGVS-token field as real
    substitutions would be misleading. Its effect is still visible as a gap
    between pid_aln (identity within the aligned core) and pid_short
    (identity over the shorter full-length sequence).
    """
    # A trailing "*" is a stop-codon marker from translation, not a residue.
    seq = seq.rstrip("*")
    seq_ref = seq_ref.rstrip("*")

    hsp = _best_hsp(seq_ref, seq)
    if hsp is None:
        return "", 0.0, 0.0

    a, b = hsp["sseq"], hsp["qseq"]  # a = ref (gapped), b = query (gapped)

    # An HSP doesn't have to start at seq_ref's first residue: ref_pos must be
    # seeded at the HSP's actual (1-based) start, or every position reported
    # from here on is silently offset.
    ref_pos = hsp["sstart"] - 1

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

    shorter = min(len(seq), len(seq_ref))
    pid_aln = round(100 * hsp["nident"] / hsp["length"], 2) if hsp["length"] > 0 else 0.0
    pid_short = round(100 * hsp["nident"] / shorter, 2) if shorter > 0 else 0.0

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
