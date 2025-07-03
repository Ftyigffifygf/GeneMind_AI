# genemind_ai/app/modules/cancer_anomaly_detector.py

def detect_anomaly(gene_sequence_str):
    """Simulates detection of abnormal gene sequences using a pre-trained ML model.
    In a real application, this would involve loading and running an actual ML model.
    """
    # Dummy logic for demonstration
    if "CANCER_MARKER_SEQ" in gene_sequence_str: # Changed marker to avoid accidental matches
        return {"is_anomalous": True, "anomaly_type": "Known Cancer Marker", "confidence": 0.95}
    elif len(gene_sequence_str) % 7 == 0 and "CANCER_MARKER_SEQ" not in gene_sequence_str: # Changed modulo and added exclusion
        return {"is_anomalous": True, "anomaly_type": "Unusual Length Pattern", "confidence": 0.70}
    else:
        return {"is_anomalous": False, "anomaly_type": "None", "confidence": 0.99}


