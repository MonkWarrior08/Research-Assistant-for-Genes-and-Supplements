import os
import streamlit as st
import openai
from Bio import Entrez
from Bio import Medline
import requests
from bs4 import BeautifulSoup
import pandas as pd
from dotenv import load_dotenv
import time
import re
import json

# Load environment variables - only if not already loaded
if not os.getenv("OPENAI_API_KEY"):
    load_dotenv()

# Set up OpenAI API key
openai.api_key = os.getenv("OPENAI_API_KEY")



# Set up NCBI email for Entrez
Entrez.email = os.getenv("ENTREZ_EMAIL")

# Configure the page - but don't override main settings
st.title("🧪 Supplement Research Assistant")
st.markdown("""
This AI assistant searches the National Center for Biotechnology Information (NCBI) 
for scientific research on supplements related to your health goals and provides a detailed analysis.
""")

# Helper text for better searches
st.markdown("""
💡 **Tip**: For best results, be specific about what you're looking for.
""")

# Functions for NCBI interaction

def search_pubmed(query, max_results=20):
    """Search PubMed for the given query and return a list of IDs."""
    # Construct a targeted search term using the original query without expansion
    search_term = f'({query}) AND (supplement OR supplementation OR nutraceutical OR herbal OR dietary OR "natural compound" OR "natural product" OR vitamin OR mineral OR nootropic OR "cognitive enhancer")'
    
    handle = Entrez.esearch(
        db="pubmed",
        term=search_term,
        retmax=max_results,
        sort="relevance",
        field="title/abstract" # Focus on searching in title and abstract for more relevant results
    )
    results = Entrez.read(handle)
    handle.close()
    return results["IdList"]

def fetch_abstracts(id_list):
    """Fetch abstracts for a list of PubMed IDs."""
    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    return list(records)

def extract_abstract_info(records, original_query):
    """Extract abstract information from PubMed records for further processing."""
    results = []
    
    # Process the original query to extract key terms
    query_terms = set(original_query.lower().split())
    
    for record in records:
        if "AB" in record:  # AB is the field for abstract
            title = record.get("TI", "No title available")
            abstract = record.get("AB", "No abstract available")
            authors = record.get("AU", ["Unknown"])
            if len(authors) > 3:
                authors = authors[:3] + ["et al."]
            authors_str = ", ".join(authors)
            journal = record.get("TA", "Unknown journal")
            year = record.get("DP", "Unknown date")[:4]  # Extract year from date
            pmid = record.get("PMID", "Unknown PMID")
            
            # We no longer extract predefined supplements
            # Instead, we'll use the actual content from the paper
            supplements = []
            
            # Score the relevance of this abstract without using expanded terms
            relevance_score = score_relevance(abstract, title, supplements, original_query, query_terms)
            
            # Check if abstract is relevant enough (minimum threshold and must contain relevant context)
            # Higher threshold (5) to ensure real relevance
            if relevance_score >= 5 and is_context_relevant(abstract, title, query_terms):
                results.append({
                    "title": title,
                    "abstract": abstract,
                    "authors": authors_str,
                    "journal": journal,
                    "year": year,
                    "pmid": pmid,
                    "supplements": supplements,
                    "relevance": relevance_score
                })
    
    # Sort results by relevance score (highest first)
    results.sort(key=lambda x: x["relevance"], reverse=True)
    return results

def is_context_relevant(abstract, title, context_terms):
    """Check if the abstract and title are relevant to the context terms."""
    combined_text = (abstract + " " + title).lower()
    
    # Must have at least 2 context terms present
    context_matches = sum(1 for term in context_terms if term in combined_text)
    
    # For stricter matching, check for key phrases
    if "improve focus" in combined_text or "cognitive function" in combined_text or "brain function" in combined_text:
        return True
        
    # Require at least 2 context matches to be considered relevant
    return context_matches >= 2

def score_relevance(abstract, title, supplements, original_query, context_terms):
    """Score the relevance of an abstract based on content and context."""
    score = 0
    combined_text = (abstract + " " + title).lower()
    
    # Removed supplements scoring since we're no longer using predefined supplements
    
    # Check for supplement mentions in title (higher weight)
    lower_title = title.lower()
    if any(term in lower_title for term in ["supplement", "vitamin", "mineral", "herb", "nutraceutical", "extract", "nootropic"]):
        score += 5
    
    # Check for exact query match (highest priority)
    if original_query.lower() in combined_text:
        score += 10
    
    # Check for context relevance
    context_matches = sum(1 for term in context_terms if term in combined_text)
    score += context_matches * 2
    
    # Check for effectiveness terms in context
    effectiveness_terms = ["effective", "efficacy", "improve", "enhance", "increase", "benefit", "significant", "result", "effect"]
    effectiveness_score = sum(1 for term in effectiveness_terms if term in combined_text)
    score += effectiveness_score
    
    # Check for research quality indicators
    quality_terms = ["randomized", "controlled trial", "double-blind", "placebo-controlled", "systematic review", "meta-analysis"]
    score += sum(2 for term in quality_terms if term in combined_text)
    
    return max(0, score)  # Ensure score doesn't go negative


def generate_analysis(abstracts, query, model="gpt4-o"):
    """Generate analysis using OpenAI API based on abstracts."""
    if not abstracts:
        return "No relevant research found."
    
    abstracts_text = ""
    for i, abstract in enumerate(abstracts):
        abstracts_text += f"[{i+1}] Title: {abstract['title']}\n"
        abstracts_text += f"    Authors: {abstract['authors']}\n"
        abstracts_text += f"    Journal: {abstract['journal']}\n"
        abstracts_text += f"    Abstract: {abstract['abstract']}\n\n"
    
    system_prompt = "You are a scientific research assistant specializing EXCLUSIVELY in nutritional supplements and their effects on health. Your task is to analyze research data and provide detailed, evidence-based analysis ONLY for supplements mentioned in the abstracts. NEVER suggest dietary changes, lifestyle modifications, or therapeutic interventions that are not supplements. Focus SOLELY on supplement-based interventions."
    
    user_prompt = f"""Based on the following research abstracts about supplements for '{query}', provide a detailed analysis FOCUSING EXCLUSIVELY ON SUPPLEMENTS.
    
    Your analysis should:
    
    1. Create a SEPARATE SECTION FOR EACH SUPPLEMENT mentioned in the abstracts
    2. For each supplement section, include:
       - Potential benefits specific to the health goal
       - Quality and strength of evidence (clinical trials, animal studies, in vitro, etc.)
       - Recommended dosages and safety considerations
       - Potential mechanisms of action explaining how the supplement works
    3. Note any conflicts between studies regarding effectiveness or recommended dosages
    
    IMPORTANT RESTRICTIONS:
    - DO NOT include any dietary recommendations or food-based interventions unless they are packaged as supplements
    - DO NOT suggest therapies, lifestyle changes, or non-supplement interventions
    - DO NOT write a conclusion or comparison section at the end
    - FOCUS ONLY on nutritional supplements, vitamins, minerals, herbs, and similar supplemental compounds
    
    Your report should be factual, evidence-based, and avoid exaggerating benefits or downplaying risks.
    If the evidence is inconclusive or limited for any supplement, clearly state this fact.
    
    RESEARCH ABSTRACTS:
    {abstracts_text}"""
    
    try:
        response = openai.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ],
            temperature=0.5,
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error generating analysis: {str(e)}"

# User input
with st.sidebar:
    st.header("Search Parameters")
    query = st.text_input("Supplement Input", placeholder="e.g. improve focus and mental clarity")
    search_btn = st.button("Research Supplements", type="primary")

if search_btn:
    if not query:
        st.error("Please enter a health goal to research.")
    elif len(query.split()) < 3:
        st.warning("Be more specific. At least 3 words to describe your health goal for better results.")
    elif not os.getenv("OPENAI_API_KEY"):
        st.error("OpenAI API key is missing. Please add it to your .env file.")
    else:
        with st.spinner("Searching for relevant research papers..."):
            # Search PubMed
            id_list = search_pubmed(query, 20)
            
            if not id_list:
                st.warning("No research papers found for the given query. Try a different search term.")
            else:
                # Fetch and process abstracts
                progress_bar = st.progress(0)
                records = fetch_abstracts(id_list)
                abstracts = extract_abstract_info(records, query)
                progress_bar.progress(50)
                
                if not abstracts:
                    st.warning("Found papers but no relevant supplement information. Try a different search term.")
                else:
                    # Generate analysis
                    st.success(f"Found {len(abstracts)} relevant research papers on supplements for '{query}'.")
                    
                    # Create tabs for Research Papers and Analysis
                    tabs = st.tabs(["Research Papers", "Analysis"])
                    
                    # Generate analysis in the background
                    with st.spinner("Generating analysis from research abstracts..."):
                        analysis = generate_analysis(abstracts, query, "gpt-4o")
                        progress_bar.progress(100)
                    
                    # Research Papers tab
                    with tabs[0]:
                        st.header(f"Research Papers for '{query}'")
                        for i, abstract in enumerate(abstracts, 1):
                            with st.expander(f"{i}. {abstract['title']}"):
                                st.markdown(f"**Authors:** {abstract['authors']}")
                                st.markdown(f"**Journal:** {abstract['journal']} ({abstract['year']})")
                                if abstract.get('supplements'):
                                    st.markdown(f"**Detected supplements:** {', '.join(abstract['supplements'])}")
                                st.markdown("### Abstract")
                                st.markdown(abstract['abstract'])
                                st.markdown(f"**PubMed ID:** {abstract['pmid']}")
                                st.markdown(f"[View on PubMed](https://pubmed.ncbi.nlm.nih.gov/{abstract['pmid']}/)")
                    
                    # Analysis tab
                    with tabs[1]:
                        st.header(f"Analysis for '{query}'")
                        st.markdown(analysis)

# Sidebar with information
with st.sidebar:
    st.subheader("About Supplement Research")
    st.markdown("""
    This assistant uses:
    - NCBI/PubMed API for scientific research
    - OpenAI gpt4-o for analysis
    
    **Note:** Always consult with a healthcare professional before starting any supplement regimen.
    """)

    
    # Add a section for API key input if not set
    if not os.getenv("OPENAI_API_KEY"):
        st.subheader("OpenAI API Key")
        api_key = st.text_input("Enter your OpenAI API key:", type="password")
        if api_key:
            os.environ["OPENAI_API_KEY"] = api_key
            st.success("API key set!")

# User interface
def main():
    # ... existing code ...
    
    # Check for API key
    if not os.getenv("OPENAI_API_KEY"):
        st.subheader("OpenAI API Key")
        api_key = st.text_input("Enter your OpenAI API key:", type="password")
        if api_key:
            os.environ["OPENAI_API_KEY"] = api_key
            st.success("API key set!")

# ... existing code ... 