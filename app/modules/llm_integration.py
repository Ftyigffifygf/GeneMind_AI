# genemind_ai/app/modules/llm_integration.py

def explain_genetic_sequence(sequence, mutations=None):
    """Simulates using an LLM to explain genetic sequences and mutations.
    In a real application, this would involve calling an external LLM API.
    """
    explanation = f"The genetic sequence provided is: {sequence}. "
    if mutations:
        explanation += "Detected mutations include: "
        for mut in mutations:
            explanation += f'{mut["type"]} at position {mut["position"]} ('
            if mut["type"] == "substitution":
                explanation += f'from {mut["reference_base"]} to {mut["patient_base"]}'
            elif mut["type"] == "insertion":
                explanation += f'inserted {mut["inserted_base"]}'
            elif mut["type"] == "deletion":
                explanation += f'deleted {mut["deleted_base"]}'
            explanation += '). '
    explanation += "This is a simplified explanation for demonstration purposes."
    return explanation


