# genemind_ai/cli.py

import argparse
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Import modules
from app.modules.gene_loader import fetch_gene_sequence
from app.modules.genome_mapping import align_sequences
from app.modules.mutation_detector import detect_mutations, correct_mutations
from app.modules.gene_creator import create_synthetic_gene
from app.modules.animal_genome_loader import load_animal_genome_dataset
from app.modules.cancer_anomaly_detector import detect_anomaly
from app.modules.llm_integration import explain_genetic_sequence
from app.modules.security_check import assess_genome_health
from app.modules.dna_converter import dna_to_gene_info
from app.modules.molecular_format_analysis import analyze_sequence

def main():
    parser = argparse.ArgumentParser(description="GeneMindAI CLI Tool")
    parser.add_argument("--api-key", type=str, default=os.getenv("API_KEY"), help="API Key for authentication")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Gene Loading
    parser_gene_loader = subparsers.add_parser("fetch-gene", help="Fetch gene sequence from NCBI")
    parser_gene_loader.add_argument("gene_id", type=str, help="NCBI Gene ID")

    # Genome Mapping
    parser_align = subparsers.add_parser("align-sequences", help="Align two DNA sequences")
    parser_align.add_argument("seq1", type=str, help="First DNA sequence")
    parser_align.add_argument("seq2", type=str, help="Second DNA sequence")

    # Mutation Detection
    parser_detect_mutations = subparsers.add_parser("detect-mutations", help="Detect mutations between two sequences")
    parser_detect_mutations.add_argument("reference_seq", type=str, help="Reference DNA sequence")
    parser_detect_mutations.add_argument("patient_seq", type=str, help="Patient DNA sequence")

    # Mutation Correction (Mock)
    parser_correct_mutations = subparsers.add_parser("correct-mutations", help="Correct mutations in a patient sequence (mock)")
    parser_correct_mutations.add_argument("patient_seq", type=str, help="Patient DNA sequence")
    parser_correct_mutations.add_argument("--mutations", type=str, help="JSON string of mutations to correct")

    # Synthetic Gene Creation
    parser_create_gene = subparsers.add_parser("create-synthetic-gene", help="Create a synthetic gene")
    parser_create_gene.add_argument("--length", type=int, default=1000, help="Length of the synthetic gene")

    # Animal Genome Loader
    parser_load_animal = subparsers.add_parser("load-animal-dataset", help="Load animal genome dataset (mock)")

    # Cancer & Anomaly Detection
    parser_detect_anomaly = subparsers.add_parser("detect-anomaly", help="Detect anomalies in a gene sequence (mock)")
    parser_detect_anomaly.add_argument("gene_sequence", type=str, help="Gene sequence to analyze")

    # LLM Integration
    parser_explain_sequence = subparsers.add_parser("explain-sequence", help="Explain genetic sequence using LLM (mock)")
    parser_explain_sequence.add_argument("sequence", type=str, help="Genetic sequence to explain")
    parser_explain_sequence.add_argument("--mutations", type=str, help="JSON string of mutations for explanation")

    # Security Check
    parser_assess_health = subparsers.add_parser("assess-health", help="Assess genome health (mock)")
    parser_assess_health.add_argument("gene_sequence", type=str, help="Gene sequence to assess")

    # DNA to Gene Info Converter
    parser_dna_to_gene_info = subparsers.add_parser("dna-to-gene-info", help="Convert DNA to gene info")
    parser_dna_to_gene_info.add_argument("dna_sequence", type=str, help="DNA sequence to convert")

    # Molecular Format Analysis
    parser_analyze_sequence = subparsers.add_parser("analyze-sequence", help="Analyze molecular format of a DNA sequence")
    parser_analyze_sequence.add_argument("dna_sequence", type=str, help="DNA sequence to analyze")

    args = parser.parse_args()

    # Simple API Key check for CLI (can be enhanced)
    if args.api_key != os.getenv("API_KEY") and args.command not in ["help"]:
        print("Error: Invalid API Key.")
        return

    if args.command == "fetch-gene":
        result = fetch_gene_sequence(args.gene_id)
        print(result)
    elif args.command == "align-sequences":
        result = align_sequences(args.seq1, args.seq2)
        print(result)
    elif args.command == "detect-mutations":
        result = detect_mutations(args.reference_seq, args.patient_seq)
        print(result)
    elif args.command == "correct-mutations":
        import json
        mutations = json.loads(args.mutations) if args.mutations else []
        result = correct_mutations(args.patient_seq, mutations)
        print(result)
    elif args.command == "create-synthetic-gene":
        result = create_synthetic_gene(args.length)
        print(result)
    elif args.command == "load-animal-dataset":
        result = load_animal_genome_dataset()
        print(result)
    elif args.command == "detect-anomaly":
        result = detect_anomaly(args.gene_sequence)
        print(result)
    elif args.command == "explain-sequence":
        import json
        mutations = json.loads(args.mutations) if args.mutations else None
        result = explain_genetic_sequence(args.sequence, mutations)
        print(result)
    elif args.command == "assess-health":
        result = assess_genome_health(args.gene_sequence)
        print(result)
    elif args.command == "dna-to-gene-info":
        result = dna_to_gene_info(args.dna_sequence)
        print(result)
    elif args.command == "analyze-sequence":
        result = analyze_sequence(args.dna_sequence)
        print(result)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()


