# genemind_ai/app/modules/molecular_format_analysis.py

from Bio.Seq import Seq
from Bio.Data import CodonTable

def analyze_sequence(dna_sequence_str):
    """Analyzes a DNA sequence for various molecular formats."""
    dna_sequence = Seq(dna_sequence_str)

    # GC Content
    if len(dna_sequence) == 0:
        gc_content = 0.0
    else:
        gc_content = (dna_sequence.count("G") + dna_sequence.count("C")) / len(dna_sequence) * 100

    # Codon Usage
    codon_usage = {}
    for i in range(0, len(dna_sequence) - len(dna_sequence) % 3, 3):
        codon = str(dna_sequence[i:i+3])
        codon_usage[codon] = codon_usage.get(codon, 0) + 1

    # Protein Translation
    try:
        # Using the standard genetic code
        protein_sequence = str(dna_sequence.translate(table=1))
    except Exception as e:
        protein_sequence = f"Translation error: {e}"

    return {
        "sequence": dna_sequence_str,
        "gc_content": gc_content,
        "codon_usage": codon_usage,
        "protein_translation": protein_sequence
    }


