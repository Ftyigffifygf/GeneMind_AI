# genemind_ai/tests/test_llm_integration.py

import unittest
from app.modules.llm_integration import explain_genetic_sequence

class TestLLMIntegration(unittest.TestCase):

    def test_explain_genetic_sequence_no_mutations(self):
        sequence = "ATGCGT"
        explanation = explain_genetic_sequence(sequence)
        self.assertIn(sequence, explanation)
        self.assertNotIn("mutations", explanation)

    def test_explain_genetic_sequence_with_mutations(self):
        sequence = "ATGCGT"
        mutations = [
            {"type": "substitution", "position": 3, "reference_base": "C", "patient_base": "G"},
            {"type": "insertion", "position": 5, "inserted_base": "A"}
        ]
        explanation = explain_genetic_sequence(sequence, mutations)
        self.assertIn(sequence, explanation)
        self.assertIn("substitution at position 3 (from C to G)", explanation)
        self.assertIn("insertion at position 5 (inserted A)", explanation)

if __name__ == "__main__":
    unittest.main()


