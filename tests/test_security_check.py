# genemind_ai/tests/test_security_check.py

import unittest
from app.modules.security_check import assess_genome_health

class TestSecurityCheck(unittest.TestCase):

    def test_assess_genome_health_cancer_marker(self):
        gene_sequence = "ATG_CANCER_MARKER_GTC"
        result = assess_genome_health(gene_sequence)
        self.assertFalse(result["is_healthy"])
        self.assertEqual(result["reason"], "Known cancer marker detected.")

    def test_assess_genome_health_pathogen_sequence(self):
        gene_sequence = "ATG_PATHOGEN_SEQUENCE_GTC"
        result = assess_genome_health(gene_sequence)
        self.assertFalse(result["is_healthy"])
        self.assertEqual(result["reason"], "Pathogen sequence detected.")

    def test_assess_genome_health_short_sequence(self):
        gene_sequence = "ATGC"
        result = assess_genome_health(gene_sequence)
        self.assertFalse(result["is_healthy"])
        self.assertEqual(result["reason"], "Abnormally short sequence.")

    def test_assess_genome_health_healthy(self):
        gene_sequence = "ATGCGT" * 10  # A sequence longer than 50
        result = assess_genome_health(gene_sequence)
        self.assertTrue(result["is_healthy"])
        self.assertEqual(result["reason"], "No obvious health risks detected (simulated).")

if __name__ == "__main__":
    unittest.main()


