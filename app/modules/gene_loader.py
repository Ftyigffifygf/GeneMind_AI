# genemind_ai/app/modules/gene_loader.py

from Bio import Entrez
from Bio import SeqIO
from app.modules.performance_utils import cached

Entrez.email = "your.email@example.com"  # Always tell NCBI who you are

@cached()
def fetch_gene_sequence(gene_id):
    """Fetches a gene sequence from the NCBI Nucleotide database."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        return {"error": str(e)}


