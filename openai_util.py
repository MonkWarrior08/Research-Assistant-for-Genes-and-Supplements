import os
from typing import Dict, List, Any, Optional
import openai
from dotenv import load_dotenv

load_dotenv()

openai.api_key = os.getenv("OPENAI_API_KEY")
if not openai.api_key:
    print("openai api key is not set")

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
            model="chatgpt-4o-latest",
            messages=[
                {"role": "system", "content": "You are a helpful scientific research assistant specializing in genetics and genomics with expertise in SNP analysis and genotype-phenotype associations."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4
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
    
    return f"""
    Based on the following research papers about {query_description}, provide a comprehensive analysis 
    focused on personal health implications and self-understanding. Do not include any title or conclusion sections.

    SECTION 1: PERSONAL GENETIC INSIGHTS
    - What the {gene_name} gene does in your body and how it affects your everyday functioning
    - How variations in this gene might influence your health, personality traits, or physical characteristics
    - Practical implications for your daily life, health management, and lifestyle choices
    - Potential risks or protective factors associated with your genetic profile
    - How environmental factors might interact with this gene to affect your wellbeing
    - Personalized context for understanding your unique genetic makeup
    
    SECTION 2: PRACTICAL APPLICATIONS & CONSIDERATIONS
    - Actionable lifestyle modifications that might be beneficial based on this genetic information
    - Nutrition, exercise, or environmental factors that may be particularly relevant
    - Questions to discuss with healthcare providers about this genetic information
    - Common misconceptions about this gene and how to properly interpret your genetic data
    - How this information fits into the broader context of your overall health
   
    {create_snp_section(snp_id, genotype)}
    
    Papers:
    {papers_text}
    
    Format your response with section headings, subheadings, and bullet points for readability.
    Use accessible, non-technical language wherever possible while maintaining accuracy.
    Focus on practical, personalized insights rather than academic discussions.
    Do not include any introduction, title, or conclusion sections in your response.
    """

def create_snp_section(snp_id: str, genotype: str) -> str:
    if not snp_id:
        return ""
    
    genotype_text = f"with genotype {genotype}" if genotype else ""
    
    return f"""
    SECTION 3: YOUR SNP PROFILE
    - Personal implications of having SNP {snp_id} {genotype_text}:
      - How this specific genetic variation might affect your body and health
      - What your genotype potentially means for your personal traits or health risks
      - How common this variation is among people with your ancestry or background
      - Lifestyle factors that might be especially important with your genetic profile
      - Practical ways to work with your genetic predispositions for optimal wellbeing
    """