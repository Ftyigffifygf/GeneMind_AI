# genemind_ai/tests/test_animal_genome_loader.py

import unittest
from app.modules.animal_genome_loader import load_animal_genome_dataset

class TestAnimalGenomeLoader(unittest.TestCase):

    def test_load_animal_genome_dataset(self):
        result = load_animal_genome_dataset()
        self.assertIn("message", result)
        self.assertIn("dataset_info", result)
        self.assertEqual(result["message"], "Animal genome dataset loaded (simulated).")
        self.assertIsInstance(result["dataset_info"], dict)
        self.assertIn("species", result["dataset_info"])
        self.assertIn("num_records", result["dataset_info"])
        self.assertIn("source", result["dataset_info"])

if __name__ == "__main__":
    unittest.main()


