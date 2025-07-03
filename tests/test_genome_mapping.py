# genemind_ai/tests/test_genome_mapping.py

import unittest
from app.modules.genome_mapping import align_sequences

class TestGenomeMapping(unittest.TestCase):

    def test_align_sequences_perfect_match(self):
        seq1 = "ATGCGT"
        seq2 = "ATGCGT"
        result = align_sequences(seq1, seq2)
        self.assertIn("score", result)
        self.assertEqual(result["score"], 6.0) # Assuming default match score of 1
        self.assertEqual(result["aligned_sequence1"], "ATGCGT")
        self.assertEqual(result["aligned_sequence2"], "ATGCGT")

    def test_align_sequences_with_mismatch(self):
        seq1 = "ATGCGT"
        seq2 = "ATGGGT"
        result = align_sequences(seq1, seq2)
        self.assertIn("score", result)
        # Score will depend on aligner parameters, but should be less than perfect
        self.assertLess(result["score"], 6.0)
        self.assertEqual(result["aligned_sequence1"], "ATGCGT")
        self.assertEqual(result["aligned_sequence2"], "ATGGGT")

    def test_align_sequences_with_gap(self):
        seq1 = "ATGCGT"
        seq2 = "ATCGT"
        result = align_sequences(seq1, seq2)
        self.assertIn("score", result)
        self.assertLess(result["score"], 6.0)
        # For global alignment, the shorter sequence will be padded with gaps
        self.assertEqual(len(result["aligned_sequence1"]), len(result["aligned_sequence2"]))
        self.assertIn("-", result["aligned_sequence1"] + result["aligned_sequence2"])

    def test_align_sequences_no_alignment(self):
        seq1 = "AAAAA"
        seq2 = "TTTTT"
        result = align_sequences(seq1, seq2)
        self.assertIn("score", result)
        self.assertLess(result["score"], 0.0) # Expecting a negative score for no match

if __name__ == "__main__":
    unittest.main()


