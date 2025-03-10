import os
from typing import Dict, List, Any, Optional
import openai
from dotenv import load_dotenv

load_dotenv()

# Set the API key for the openai package (v0.28.1)
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    print("openai api key is not set")
else:
    openai.api_key = api_key

def analyze_papers(papers: List[Dict[str, Any]], gene_name: str, analysis_type: str = "comprehensive", snp_id: str = "", genotype: str = "") -> str:
    if not papers:
        return "No paper to analyze."
    
    try:
        paper_texts = []
        for i, paper in enumerate(papers[:20]):
            paper_text = f"PAPER {i+1}:\n"
            paper_text += f"TITLE: {paper['title']}\n"
            paper_text += f"AUTHORS: {paper['author']}\n"
            paper_text += f"JOURNAL: {paper['journal']} ({paper['publication_date']})\n"
            paper_text += f"ABSTRACT: {paper['abstract']}\n\n"
            paper_texts.append(paper_text)
        
        all_paper_text = "\n".join(paper_texts)

        prompt = create_comprehensive_prompt(gene_name, all_paper_text, snp_id, genotype)

        
        response = openai.chat.completions.create(
            model="gpt-4o",  
            messages=[
                {"role": "system", "content": "You are a helpful scientific research assistant specializing in genetics and genomics with expertise in SNP analysis and genotype-phenotype associations."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.5
        )

        return response.choices[0].message.content
    
    except Exception as e:
        print(f"Error during openai analysis: {e}")
        return f"Error performing analysis: {str(e)}"
    
def create_comprehensive_prompt(gene_name: str, papers_text: str, snp_id: str = "", genotype: str = "") -> str:
    query_description = f"gene {gene_name}"
    if snp_id:
        query_description = f"gene {gene_name} and SNP {snp_id}"
        if genotype:
            query_description += f" with genotype {genotype}"
    
    # Create SNP-specific context if applicable
    snp_context = ""
    if snp_id and genotype:
        snp_context = f" and your specific {snp_id} {genotype} genotype (both alleles)"
    elif snp_id:
        snp_context = f" and your SNP {snp_id}"
    
    return f"""
    Based on the following research papers about {query_description}, provide a comprehensive detailed scientific analysis 
    focused on personal health implications and self-understanding. Do not include any title or conclusion sections.

    SECTION 1: PERSONAL GENETIC INSIGHTS
    - What the {gene_name} gene does in your body and how it affects your everyday functioning{snp_context}
    - How your specific genetic variation{snp_context} might influence your health, personality traits, or physical characteristics
  
    SECTION 2: PRACTICAL APPLICATIONS & CONSIDERATIONS
    - First highlight any practical applications or recommendations specifically mentioned in the research papers for your genotype{snp_context}
    - Focus on domain-specific factors directly related to the {gene_name} gene's function and your specific genetic variation(e.g., if it relates to social behavior, focus on social interactions) 
    
    Papers:
    {papers_text}
    
    Format your response with section headings, subheadings, and bullet points for readability.
    Use accessible and technical language wherever possible while maintaining accuracy.
    Focus on practical, personalized insights with deep scientific analysis.
    Do not include any introduction, title, or conclusion sections in your response.
    """