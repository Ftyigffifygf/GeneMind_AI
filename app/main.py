# genemind_ai/app/main.py

from fastapi import FastAPI, Depends, HTTPException, status
from fastapi.security import APIKeyHeader
from dotenv import load_dotenv
import os

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

app = FastAPI(title="GeneMindAI API", description="Advanced Genomic Analysis Tool")

API_KEY = os.getenv("API_KEY")
API_KEY_NAME = "X-API-Key"
api_key_header = APIKeyHeader(name=API_KEY_NAME, auto_error=True)

async def get_api_key(api_key: str = Depends(api_key_header)):
    if api_key == API_KEY:
        return api_key
    raise HTTPException(
        status_code=status.HTTP_403_FORBIDDEN,
        detail="Could not validate credentials",
    )

@app.get("/", summary="Root endpoint", tags=["General"])
async def read_root():
    return {"message": "Welcome to GeneMindAI API!"}

# Endpoints for each module

@app.get("/gene_loader/fetch_sequence/{gene_id}", summary="Fetch Gene Sequence", tags=["Gene Loading"])
async def fetch_sequence_endpoint(gene_id: str, api_key: str = Depends(get_api_key)):
    return fetch_gene_sequence(gene_id)

@app.post("/genome_mapping/align_sequences", summary="Align DNA Sequences", tags=["Genome Mapping"])
async def align_sequences_endpoint(seq1: str, seq2: str, api_key: str = Depends(get_api_key)):
    return align_sequences(seq1, seq2)

@app.post("/mutation_detection/detect_mutations", summary="Detect Mutations", tags=["Mutation Detection and Correction"])
async def detect_mutations_endpoint(reference_seq: str, patient_seq: str, api_key: str = Depends(get_api_key)):
    return detect_mutations(reference_seq, patient_seq)

@app.post("/mutation_detection/correct_mutations", summary="Correct Mutations (Mock)", tags=["Mutation Detection and Correction"])
async def correct_mutations_endpoint(patient_seq: str, mutations: list, api_key: str = Depends(get_api_key)):
    return correct_mutations(patient_seq, mutations)

@app.get("/gene_creator/create_synthetic_gene", summary="Create Synthetic Gene (Mock)", tags=["Synthetic Gene Creation"])
async def create_synthetic_gene_endpoint(length: int = 1000, api_key: str = Depends(get_api_key)):
    return {"synthetic_gene": create_synthetic_gene(length)}

@app.get("/animal_genome_loader/load_dataset", summary="Load Animal Genome Dataset (Mock)", tags=["Animal Genome Dataset Loader"])
async def load_animal_genome_dataset_endpoint(api_key: str = Depends(get_api_key)):
    return load_animal_genome_dataset()

@app.post("/cancer_anomaly_detection/detect_anomaly", summary="Detect Cancer/Anomaly (Mock)", tags=["Cancer & Anomaly Detection"])
async def detect_anomaly_endpoint(gene_sequence: str, api_key: str = Depends(get_api_key)):
    return detect_anomaly(gene_sequence)

@app.post("/llm_integration/explain_sequence", summary="Explain Genetic Sequence (LLM Mock)", tags=["LLM API Integration"])
async def explain_sequence_endpoint(sequence: str, mutations: list = None, api_key: str = Depends(get_api_key)):
    return {"explanation": explain_genetic_sequence(sequence, mutations)}

@app.post("/security_check/assess_health", summary="Assess Genome Health (Mock)", tags=["Security Check"])
async def assess_health_endpoint(gene_sequence: str, api_key: str = Depends(get_api_key)):
    return assess_genome_health(gene_sequence)

@app.post("/dna_converter/dna_to_gene_info", summary="DNA to Gene Info Converter", tags=["DNA to Gene Info Converter"])
async def dna_to_gene_info_endpoint(dna_sequence: str, api_key: str = Depends(get_api_key)):
    return dna_to_gene_info(dna_sequence)

@app.post("/molecular_format_analysis/analyze_sequence", summary="Analyze Molecular Format", tags=["Molecular Format Analysis"])
async def analyze_sequence_endpoint(dna_sequence: str, api_key: str = Depends(get_api_key)):
    return analyze_sequence(dna_sequence)


