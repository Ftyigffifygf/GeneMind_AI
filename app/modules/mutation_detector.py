# genemind_ai/app/modules/mutation_detector.py

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def detect_mutations(reference_seq_str, patient_seq_str):
    """Detects mutations by comparing a patient gene sequence with a reference sequence.
    Currently focuses on substitutions for simplicity and robustness.
    """
    mutations = []
    min_len = min(len(reference_seq_str), len(patient_seq_str))

    for i in range(min_len):
        if reference_seq_str[i] != patient_seq_str[i]:
            mutations.append({
                "type": "substitution",
                "position": i,
                "reference_base": reference_seq_str[i],
                "patient_base": patient_seq_str[i]
            })
            
    # Handle simple insertions/deletions at the end if one sequence is longer
    if len(reference_seq_str) < len(patient_seq_str):
        for i in range(len(reference_seq_str), len(patient_seq_str)):
            mutations.append({"type": "insertion", "position": i, "inserted_base": patient_seq_str[i]})
    elif len(patient_seq_str) < len(reference_seq_str):
        for i in range(len(patient_seq_str), len(reference_seq_str)):
            mutations.append({"type": "deletion", "position": i, "deleted_base": reference_seq_str[i]})

    return mutations

def correct_mutations(patient_seq_str, mutations):
    """Simulates correcting mutations in a patient's gene sequence.
    This mock correction primarily handles substitutions.
    """
    corrected_seq_list = list(patient_seq_str)
    
    # Sort mutations by position in descending order to avoid index issues during correction
    mutations.sort(key=lambda x: x["position"], reverse=True)

    for mutation in mutations:
        if mutation["type"] == "substitution":
            pos = mutation["position"]
            if pos < len(corrected_seq_list) and corrected_seq_list[pos] == mutation["patient_base"]:
                corrected_seq_list[pos] = mutation["reference_base"]
        # For insertions and deletions, a simple correction is not robust without re-alignment.
        # These are left as identified but not corrected in this simplified mock.

    return "".join(corrected_seq_list)


