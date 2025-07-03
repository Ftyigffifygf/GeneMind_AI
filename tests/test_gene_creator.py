# genemind_ai/tests/test_gene_creator.py

import unittest
from app.modules.gene_creator import create_synthetic_gene

class TestGeneCreator(unittest.TestCase):

    def test_create_synthetic_gene_default_length(self):
        gene = create_synthetic_gene()
        self.assertEqual(len(gene), 1000)
        self.assertTrue(all(base in "ATCG" for base in gene))

    def test_create_synthetic_gene_custom_length(self):
        gene = create_synthetic_gene(length=500)
        self.assertEqual(len(gene), 500)
        self.assertTrue(all(base in "ATCG" for base in gene))

    def test_create_synthetic_gene_zero_length(self):
        gene = create_synthetic_gene(length=0)
        self.assertEqual(len(gene), 0)
        self.assertEqual(gene, "")

if __name__ == "__main__":
    unittest.main()


