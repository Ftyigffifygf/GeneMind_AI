# genemind_ai/app/modules/security_check.py

def assess_genome_health(gene_sequence_str):
    """Simulates scanning a given genome to assess if it is healthy/safe.
    This is a mock function for educational/research purposes.
    """
    # Dummy logic for demonstration
    if "CANCER_MARKER" in gene_sequence_str:
        return {"is_healthy": False, "reason": "Known cancer marker detected."}
    elif "PATHOGEN_SEQUENCE" in gene_sequence_str:
        return {"is_healthy": False, "reason": "Pathogen sequence detected."}
    elif len(gene_sequence_str) < 50:
        return {"is_healthy": False, "reason": "Abnormally short sequence."}
    else:
        return {"is_healthy": True, "reason": "No obvious health risks detected (simulated)."}


