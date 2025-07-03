# genemind_ai/app/modules/dna_converter.py

from Bio.Seq import Seq
from Bio.Data import CodonTable

def dna_to_gene_info(dna_sequence_str):
    """Translates DNA into mRNA, amino acids, and provides predicted protein functions (mock)."""
    dna_sequence = Seq(dna_sequence_str)

    # DNA to mRNA
    mRNA_sequence = dna_sequence.transcribe()

    # mRNA to Amino Acids (Protein Translation)
    # Using the standard genetic code (table=1)
    protein_sequence = mRNA_sequence.translate(table=1, cds=False)

    # Predicted Protein Functions (Mock/Simplified)
    # In a real application, this would involve complex bioinformatics analysis
    predicted_functions = []
    if "ATGC" in dna_sequence_str:
        predicted_functions.append("Potential regulatory region")
    if len(protein_sequence) > 50:
        predicted_functions.append("Likely a functional protein")
    else:
        predicted_functions.append("Short peptide or non-coding region")

    return {
        "dna_sequence": dna_sequence_str,
        "mRNA_sequence": str(mRNA_sequence),
        "protein_sequence": str(protein_sequence),
        "predicted_functions": predicted_functions
    }


