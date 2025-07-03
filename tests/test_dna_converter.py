# genemind_ai/tests/test_dna_converter.py

import unittest
from app.modules.dna_converter import dna_to_gene_info

class TestDnaConverter(unittest.TestCase):

    def test_dna_to_gene_info_transcription(self):
        dna_sequence = "ATGCGT"
        result = dna_to_gene_info(dna_sequence)
        self.assertEqual(result["mRNA_sequence"], "AUGCGU")

    def test_dna_to_gene_info_translation(self):
        dna_sequence = "ATGCGT"
        result = dna_to_gene_info(dna_sequence)
        self.assertEqual(result["protein_sequence"], "MR")

    def test_dna_to_gene_info_predicted_functions(self):
        dna_sequence = "ATGCGT" * 30 # Long enough sequence to produce protein > 50 amino acids
        result = dna_to_gene_info(dna_sequence)
        self.assertIn("Likely a functional protein", result["predicted_functions"])

        short_dna_sequence = "ATGC"
        short_result = dna_to_gene_info(short_dna_sequence)
        self.assertIn("Short peptide or non-coding region", short_result["predicted_functions"])

if __name__ == "__main__":
    unittest.main()


