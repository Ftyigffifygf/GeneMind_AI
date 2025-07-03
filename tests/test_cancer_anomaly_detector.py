# genemind_ai/tests/test_cancer_anomaly_detector.py

import unittest
from app.modules.cancer_anomaly_detector import detect_anomaly

class TestCancerAnomalyDetector(unittest.TestCase):

    def test_detect_anomaly_known_marker(self):
        gene_sequence = "ATG_CANCER_MARKER_SEQ_GTC"
        result = detect_anomaly(gene_sequence)
        self.assertTrue(result["is_anomalous"])
        self.assertEqual(result["anomaly_type"], "Known Cancer Marker")
        self.assertEqual(result["confidence"], 0.95)

    def test_detect_anomaly_unusual_length(self):
        gene_sequence = "ATGCATG" # Length is 7
        result = detect_anomaly(gene_sequence)
        self.assertTrue(result["is_anomalous"])
        self.assertEqual(result["anomaly_type"], "Unusual Length Pattern")
        self.assertEqual(result["confidence"], 0.70)

    def test_detect_anomaly_no_anomaly(self):
        gene_sequence = "ATGCATGC" # Length is 8
        result = detect_anomaly(gene_sequence)
        self.assertFalse(result["is_anomalous"])
        self.assertEqual(result["anomaly_type"], "None")
        self.assertEqual(result["confidence"], 0.99)

if __name__ == "__main__":
    unittest.main()


