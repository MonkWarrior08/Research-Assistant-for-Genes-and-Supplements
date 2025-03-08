import os
from typing import Dict, List, Optional, Any
from Bio import Entrez
import time
import re
import requests
from xml.etree import ElementTree
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")
if Entrez.email:
    print("NCBI_EMAIL set")
else:
    print("NCBI_EMAIL is not set")

ncbi_api_key = os.getenv("NCBI_API_KEY")
if ncbi_api_key:
    Entrez.api_key = ncbi_api_key
    print("NCBI_API_KEY set")
else:
    print("NCBI_API_KEY is not set")

def parse_gene_input(gene_input: str) -> Dict[str, str]:
    parts = gene_input.strip().split()
    result = {
        "gene": "",
        "rs_id": "",
        "genotype": "",
        "original_input": gene_input
    }
    if len(parts)>= 1:
        result["gene"] = parts[0]
    if len(parts)>= 2:
        if re.match(r'rs\d+', parts[1]):
            result["rs_id"] = parts[1]
    if len(parts) >= 3:
        if re.match(r'[ACgT]{2}', parts[2]):
            result["genotype"] = parts[2]
    return result

def search_gene_paper(gene_name: str, max_results: int=20) -> list[str]:
    gene_data = parse_gene_input(gene_name)
    search_terms = []

    if gene_data["gene"]:
        search_terms.append(f"{gene_data['gene']}[Gene/Protein]")
    if gene_data["rs_id"]:
        search_terms.append(f"{gene_data['rs_id']}[All Fileds]")
    if gene_data["genotype"]:
        allele_search = " OR ".join([f"{allele}[All Fields]" for allele in gene_data["genotype"]])
        search_terms.append(f"({allele_search})")
    
    search_term = " AND ".join(search_terms)

    try:
        search_handle = Entrez.esearch(
            db="pubmed",
            term=search_term,
            retmax=max_results,
            sort="relevance"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        return search_results["IdList"]
    except Exception as e:
        print(f"Error searching: {e}")
        return []

def fetch_paper_details(pmid_list: List[str]) -> List[Dict[str, Any]]:
    if not pmid_list:
        return []
    papers = []

    try:
        handle = Entrez.esummary(db="pubmed", id=",".join(pmid_list))
        summary_result = Entrez.read(handle)
        handle.close()

        for i, pmid in enumerate(pmid_list):
            if i >= len(summary_result):
                break

            paper_info = summary_result[i]

            paper = {
                "pmid": pmid,
                "title": paper_info.get("Title", ""),
            }