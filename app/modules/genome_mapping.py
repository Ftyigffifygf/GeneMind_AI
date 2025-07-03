# genemind_ai/app/modules/genome_mapping.py

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from app.modules.performance_utils import cached

@cached()
def align_sequences(seq1_str, seq2_str):
    """Aligns two DNA sequences and returns the alignment score and aligned sequences."""
    seq1 = Seq(seq1_str)
    seq2 = Seq(seq2_str)
    aligner = PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.mode = 'global'
    alignments = aligner.align(seq1, seq2)
    
    if alignments:
        top_alignment = alignments[0]
        return {
            "score": top_alignment.score,
            "aligned_sequence1": str(top_alignment[0]),
            "aligned_sequence2": str(top_alignment[1])
        }
    else:
        return {"error": "No alignment found."}


