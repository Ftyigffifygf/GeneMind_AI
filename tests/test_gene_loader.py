# genemind_ai/tests/test_gene_loader.py

import unittest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from app.modules.gene_loader import fetch_gene_sequence

class TestGeneLoader(unittest.TestCase):

    @patch("Bio.Entrez.efetch")
    @patch("Bio.SeqIO.read")
    def test_fetch_gene_sequence_success(self, mock_seq_read, mock_efetch):
        # Mock the Entrez.efetch and SeqIO.read calls
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_seq_read.return_value = SeqRecord(Seq("ATGC"), id="NM_000001", name="test_gene")

        gene_id = "NM_000001"
        result = fetch_gene_sequence(gene_id)

        mock_efetch.assert_called_once_with(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        mock_seq_read.assert_called_once_with(mock_handle, "fasta")
        self.assertIsInstance(result, SeqRecord)
        self.assertEqual(str(result.seq), "ATGC")
        self.assertEqual(result.id, "NM_000001")

    @patch("Bio.Entrez.efetch", side_effect=Exception("Network error"))
    def test_fetch_gene_sequence_error(self, mock_efetch):
        gene_id = "NM_000002"
        result = fetch_gene_sequence(gene_id)

        mock_efetch.assert_called_once_with(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        self.assertIn("error", result)
        self.assertEqual(result["error"], "Network error")

if __name__ == "__main__":
    unittest.main()


