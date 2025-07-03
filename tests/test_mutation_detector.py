# genemind_ai/tests/test_mutation_detector.py

import unittest
from app.modules.mutation_detector import detect_mutations, correct_mutations

class TestMutationDetector(unittest.TestCase):

    def test_detect_mutations_substitution(self):
        ref_seq = "ATGCGT"
        pat_seq = "ATGGGT"
        mutations = detect_mutations(ref_seq, pat_seq)
        expected_mutations = [{
            "type": "substitution",
            "position": 3,
            "reference_base": "C",
            "patient_base": "G"
        }]
        self.assertEqual(mutations, expected_mutations)

    def test_detect_mutations_insertion_at_end(self):
        ref_seq = "ATGC"
        pat_seq = "ATGCGT"
        mutations = detect_mutations(ref_seq, pat_seq)
        expected_mutations = [
            {"type": "insertion", "position": 4, "inserted_base": "G"},
            {"type": "insertion", "position": 5, "inserted_base": "T"}
        ]
        self.assertEqual(mutations, expected_mutations)

    def test_detect_mutations_deletion_at_end(self):
        ref_seq = "ATGCGT"
        pat_seq = "ATGC"
        mutations = detect_mutations(ref_seq, pat_seq)
        expected_mutations = [
            {"type": "deletion", "position": 4, "deleted_base": "G"},
            {"type": "deletion", "position": 5, "deleted_base": "T"}
        ]
        self.assertEqual(mutations, expected_mutations)

    def test_detect_mutations_no_mutation(self):
        ref_seq = "ATGCGT"
        pat_seq = "ATGCGT"
        mutations = detect_mutations(ref_seq, pat_seq)
        self.assertEqual(mutations, [])

    def test_correct_mutations_substitution(self):
        patient_seq = "ATGGGT"
        mutations = [{
            "type": "substitution",
            "position": 3,
            "reference_base": "C",
            "patient_base": "G"
        }]
        corrected_seq = correct_mutations(patient_seq, mutations)
        self.assertEqual(corrected_seq, "ATGCGT")

    def test_correct_mutations_no_correction(self):
        patient_seq = "ATGCGT"
        mutations = []
        corrected_seq = correct_mutations(patient_seq, mutations)
        self.assertEqual(corrected_seq, "ATGCGT")

if __name__ == "__main__":
    unittest.main()


