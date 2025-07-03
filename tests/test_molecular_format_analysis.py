# genemind_ai/tests/test_molecular_format_analysis.py

import unittest
from app.modules.molecular_format_analysis import analyze_sequence

class TestMolecularFormatAnalysis(unittest.TestCase):

    def test_analyze_sequence_gc_content(self):
        dna_sequence = "ATGCATGC"
        result = analyze_sequence(dna_sequence)
        self.assertIn("gc_content", result)
        self.assertEqual(result["gc_content"], 50.0)

    def test_analyze_sequence_codon_usage(self):
        dna_sequence = "ATGCATGCATGCATGCATGC"
        result = analyze_sequence(dna_sequence)
        self.assertIn("codon_usage", result)
        # Corrected expected codon usage based on the sequence
        self.assertEqual(result["codon_usage"], {"ATG": 2, "TGC": 1, "GCA": 1, "CAT": 2})

    def test_analyze_sequence_protein_translation(self):
        dna_sequence = "ATGCATGCATGCATGCATGC"
        result = analyze_sequence(dna_sequence)
        self.assertIn("protein_translation", result)
        # Corrected expected protein translation
        self.assertEqual(result["protein_translation"], "MHACMH")

    def test_analyze_sequence_empty(self):
        dna_sequence = ""
        result = analyze_sequence(dna_sequence)
        self.assertIn("gc_content", result)
        self.assertEqual(result["gc_content"], 0.0)
        self.assertEqual(result["codon_usage"], {})
        self.assertEqual(result["protein_translation"], "") # Biopython returns empty string for empty sequence

if __name__ == "__main__":
    unittest.main()


