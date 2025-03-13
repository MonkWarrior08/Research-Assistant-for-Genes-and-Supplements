import os
from typing import Dict, List, Any, Optional
import openai
from dotenv import load_dotenv

load_dotenv()

# Set the API key for the openai package (v0.28.1)cs
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    print("openai api key is not set")
else:
    openai.api_key = api_key

def analyze_papers(papers: List[Dict[str, Any]], gene_name: str, analysis_type: str = "comprehensive", snp_id: str = "", genotype: str = "") -> str:
    if not papers:
        return "No paper to analyze."
    
    try:
        # Process each abstract individually
        results = []
        for i, paper in enumerate(papers[:20]):
            if paper['abstract']:
                # Analyze each abstract with detailed bullet points
                abstract_analysis = analyze_single_abstract(
                    paper['abstract'], 
                    gene_name, 
                    snp_id, 
                    genotype
                )
                
                if abstract_analysis and abstract_analysis != "No relevant information found.":
                    paper_result = f"## PAPER {i+1}:\n"
                    paper_result += f"**TITLE**: {paper['title']}\n"
                    paper_result += f"**AUTHORS**: {paper['author']}\n"
                    paper_result += f"**JOURNAL**: {paper['journal']} ({paper['publication_date']})\n"
                    paper_result += f"**URL**: {paper.get('url', '')}\n\n"
                    paper_result += f"{abstract_analysis}\n\n"
                    
                    results.append(paper_result)
        
        if not results:
            return f"No relevant information about {gene_name} gene, SNP {snp_id}, or genotype {genotype} found in any of the papers."
        
        return "\n".join(results)
    
    except Exception as e:
        print(f"Error during abstract analysis: {e}")
        return f"Error analyzing abstracts: {str(e)}"

def analyze_single_abstract(abstract: str, gene_name: str, snp_id: str = "", genotype: str = "") -> str:
    """Analyze a single abstract with detailed bullet points for SNP ID and genotype information."""
    try:
        # Create a focused prompt for analyzing the abstract with bullet points
        analysis_prompt = f"""
        Analyze this scientific abstract to extract detailed information about the gene {gene_name.upper()},
        focusing specifically on the SNP {snp_id} and genotype {genotype} if mentioned.
        
        For each section (SNP and genotype), provide:
        1. Main bullet points of all relevant information found in the abstract
        2. For each main bullet point, add a sub-bullet point that explains its meaning and provides clarification
        
        Format your response as follows:
        
        ### SNP {snp_id} Information:
        • [Key finding or information about the SNP from the abstract]
          • [Explanation/clarification of what this means in simple terms]
        • [Another key finding about the SNP]
          • [Explanation/clarification]
        
        ### Genotype {genotype} Information:
        • [Key finding or information about this specific genotype]
          • [Explanation/clarification of what this means in simple terms]
        • [Another key finding about this genotype]
          • [Explanation/clarification]
        
        Important guidelines:
        - If no information is found about the SNP, include only the gene information.
        - If no information is found about the genotype, omit that section entirely.
        - If neither SNP nor genotype information is found, respond with "No relevant information found."
        - Use bullet points (•) for main points and sub-bullets for explanations.
        - Keep the language precise yet accessible.
        
        Abstract:
        {abstract}
        """
        
        response = openai.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are a scientific research assistant specializing in genetics. You extract and explain genetic information from research papers in a clear, organized way using bullet points."},
                {"role": "user", "content": analysis_prompt}
            ],
            temperature=0.1,
            max_tokens=1000
        )
        
        return response.choices[0].message.content.strip()
    
    except Exception as e:
        print(f"Error during abstract analysis: {e}")
        return f"Error analyzing abstract: {str(e)}"