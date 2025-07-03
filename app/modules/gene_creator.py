# genemind_ai/app/modules/gene_creator.py

import random

def create_synthetic_gene(length=1000):
    """Generates a random synthetic DNA sequence of a given length."""
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) for _ in range(length))


