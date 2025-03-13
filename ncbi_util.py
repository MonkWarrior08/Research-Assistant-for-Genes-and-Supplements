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

def search_gene_paper(gene_name: str, max_results: int=20, custom_date_range: Optional[str]=None) -> list[str]:
    gene_data = parse_gene_input(gene_name)
    search_terms = []
    
    # Build search terms with a balanced approach - more focused than the previous version
    # but still more inclusive than the original
    
    # Add gene search term with a broader field coverage
    if gene_data["gene"]:
        search_terms.append(f"({gene_data['gene']}[Gene Symbol] OR {gene_data['gene']}[Title/Abstract])")
    
    # Add SNP search term
    if gene_data["rs_id"]:
        search_terms.append(f"{gene_data['rs_id']}[All Fields]")
    
    # Add genotype search term
    if gene_data["genotype"]:
        allele_search = " OR ".join([f"{allele}[All Fields]" for allele in gene_data["genotype"]])
        search_terms.append(f"({allele_search})")
    
    # Combine all search terms with AND for more relevant results
    search_term = " AND ".join(search_terms)

    # Add date range to search if provided
    search_params = {
        "db": "pubmed",
        "term": search_term,
        "retmax": max_results,
        "sort": "relevance"
    }

    if custom_date_range:
        search_params["datetype"] = "pdat"  # Publication date
        search_params["mindate"], search_params["maxdate"] = custom_date_range.split(":")

    try:
        search_handle = Entrez.esearch(**search_params)
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
                "author": ", ".join(paper_info.get("AuthorList", [])),
                "journal": paper_info.get("FullJournalName", ""),
                "publication_date": paper_info.get("PubDate", ""),
                "abstract": "",
                "doi": "",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }
            papers.append(paper)
        
        for paper in papers:
            fetch_abstract(paper)
            time.sleep(0.5)
        return papers
    
    except Exception as e:
        print(f"Error fetching paper details: {e}")
        return papers
    
def fetch_abstract(paper: Dict[str, Any]) -> None:
    try:
        handle = Entrez.efetch(
            db="pubmed",
            id=paper["pmid"],
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()

        article = record["PubmedArticle"][0]["MedlineCitation"]["Article"]

        if "Abstract" in article:
            abstract_parts = article["Abstract"]["AbstractText"]

            if isinstance(abstract_parts, list):
                full_abstract = ""
                for part in abstract_parts:
                    if hasattr(part, "attributes") and "Label" in part.attributes:
                        full_abstract += f"{part.attributes['Label']}: {part}\n"
                    else:
                        full_abstract += f"{part}\n"
                paper["abstract"] = full_abstract.strip()
            else:
                paper["abstract"] = abstract_parts

        article_id = record["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
        for id_item in article_id:
            if id_item.attributes.get("IdType") == "doi":
                paper["doi"] = str(id_item)
                break
    except Exception as e:
        print(f"Error fetching abstract for {paper['pmid']}: {e}")

def get_gene_info(gene_name: str) -> Dict[str, Any]:
    gene_data = parse_gene_input(gene_name)
    gene_symbol = gene_data["gene"]
    rs_id = gene_data["rs_id"]

    results = {
        "gene_info": {},
        "snp_info": {},
        "original_query": gene_name
    }
    
    try:
        if gene_symbol:
            search_handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Name] OR {gene_symbol}[Gene Symbol]")
            search_result = Entrez.read(search_handle)
            search_handle.close()

            if search_result["IdList"]:
                gene_id = search_result["IdList"][0]
                fetch_handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                gene_record = Entrez.read(fetch_handle)
                fetch_handle.close()

                if gene_record:
                    gene_data = gene_record[0]

                    results["gene_info"]= {
                        "gene_id": gene_id,
                        "name": gene_data.get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_locus", ""),
                        "symbol": gene_data.get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_locus-tag", ""),
                        "description": gene_data.get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_desc", ""),
                        "organism": gene_data.get("Entrezgene_source", {}).get("BioSource", {}).get("BioSource_org", {}).get("Org-ref", {}).get("Org-ref_taxname", ""),
                        "aliases": gene_data.get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_syn", []),
                    }
            else:
                results["gene_info"] = {"error": f"Gene {gene_symbol} not foud"}

    except Exception as e:
        results["gene_info"] = {"error": f"Error retrieving gene info: {str(e)}"}
    
    if rs_id:
        results["snp_info"] = get_snp_info_direct(rs_id)
        genotype = gene_data.get("genotype", "")
        if genotype:
            results["snp_info"]["genotype"] = genotype
    return results

def get_snp_info_direct(rs_id: str) -> Dict[str, Any]:
    if not rs_id:
        return {"error": "No SNP Id provided"}
    
    rs_number = rs_id.lower().replace("rs", "")

    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id={rs_number}"
        if ncbi_api_key:
            url += f"&api_key={ncbi_api_key}"
        
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        root = ElementTree.fromstring(response.content)

        snp_info = {
            "rs_id": rs_id,
            "variant_type": "",
            "gene_names": [],
            "clinical_significance": "",
            "chromosome": "",
            "position": "",
            "minor_allele": "",
            "minor_allele_freq": "",
            "assembly": "",
            "alleles": []
        }

        for item in root.findall(".//Item"):
            name = item.get("Name")
            if name == "CHRPOS" and item.text:
                chrpos = item.text.split(":")
                if len(chrpos) >= 2:
                    snp_info["chromosome"] = chrpos[0]
                    snp_info["position"] = chrpos[1]
            elif name == "GENES" and len(item) > 0:
                for gene in item.findall(".//Item[@Name='NAME']"):
                    if gene.text:
                        snp_info["gene_names"].append(gene.text)
            elif name == "SNP_CLASS" and item.text:
                snp_info["variant_type"] = item.text
            elif name == "CLINICAL_SIGNIFICANCE" and item.text:
                snp_info["clinical_significance"] = item.text
            elif name == "GLOBAL_MAFS" and len(item) > 0:
                for maf_item in item.findall(".//Item[@Name='STUDY']"):
                    study_name = maf_item.text
                    if study_name and "1000Genomes" in study_name:
                        for sibling in maf_item.getparent():
                            if sibling.get("Name") == "FREQ":
                                freqs = sibling.text.split("/") if sibling.text else []
                                if len(freqs) >= 2:
                                    snp_info["minor_allele"] = freqs[0]
                                    snp_info["minor_allele_freq"] = freqs[1]
                                break
        if not snp_info["gene_names"]:
            snp_info["gene_names"] = ["Not specified"]
        if not snp_info["variant_type"]:
            snp_info["variant_type"] = "unknown"
        if not snp_info["clinical_significance"]:
            snp_info["clinical_significance"] = "Not specified"
        if not snp_info["chromosome"]:
            snp_info["chromosome"] = "Not specified"
        if not snp_info["position"]:
            snp_info["position"] = "Not specified"
        
        return snp_info
    
    except requests.exceptions.RequestException as e:
        print(f"HTTP request error for SNP {rs_id}: {e}")
        return {"error": f"Error connecting to NCBI database: {str(e)}"}
    except ElementTree.ParseError as e:
        print(f"XML parsing error for SNP {rs_id}: {e}")
        return {"error": f"Error parsing SNP data: {str(e)}"}
    except Exception as e:
        print(f"Unexpected error retrieving SNP {rs_id}: {e}")
        return {"error": f"Unexpected error: {str(e)}"}