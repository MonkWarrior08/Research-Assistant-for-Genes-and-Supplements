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
from cache_utils import cached_function

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
@cached_function
def search_pubmed(query, max_results=20):
    """Search PubMed for the given query and return a list of IDs."""
    # Expand query with common synonyms and related terms
    expanded_terms = expand_query_with_synonyms(query)
    
    # Construct a more targeted search term
    search_term = f'({expanded_terms}) AND (supplement OR supplementation OR nutraceutical OR herbal OR dietary OR "natural compound" OR "natural product" OR vitamin OR mineral OR nootropic OR "cognitive enhancer")'
    
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

def expand_query_with_synonyms(query):
    """Expand a query with synonyms and related terms for more comprehensive search."""
    # Dictionary of common health goals and their synonyms/related terms
    synonyms_dict = {
        "focus": 'focus OR attention OR concentration OR cognition OR cognitive OR "mental clarity" OR alertness OR "mental performance" OR "brain function" OR "cognitive performance" OR "cognitive enhancement" OR "mental focus" OR "attention span" OR "cognitive ability"',
        "energy": 'energy OR fatigue OR tiredness OR "mental energy" OR vitality OR stamina OR endurance OR "physical performance" OR "mental fatigue" OR "physical fatigue" OR "energy boost" OR "energy level"',
        "sleep": 'sleep OR insomnia OR "sleep quality" OR "sleep duration" OR "sleep disorder" OR "sleep disturbance" OR "sleep onset" OR "sleep maintenance" OR "deep sleep" OR "REM sleep" OR "sleep efficiency" OR "sleep architecture"',
        "stress": 'stress OR anxiety OR tension OR "stress relief" OR "stress reduction" OR relaxation OR "cortisol level" OR "stress hormone" OR "stress management" OR "anxiety reduction" OR calm OR "mental stress" OR "chronic stress"',
        "inflammation": 'inflammation OR "inflammatory marker" OR "anti-inflammatory" OR "chronic inflammation" OR cytokine OR "inflammatory response" OR "systemic inflammation" OR "inflammatory cytokine" OR "inflammation reduction"',
        "memory": 'memory OR "memory enhancement" OR "memory formation" OR "memory recall" OR "memory retention" OR "short-term memory" OR "long-term memory" OR "working memory" OR "memory function" OR "memory improvement" OR "memory consolidation"',
        "mood": 'mood OR depression OR "mood enhancement" OR "mood improvement" OR "mood regulation" OR "emotional wellbeing" OR "positive mood" OR "mood stabilization" OR "mood disorder" OR "mood support" OR happiness OR "emotional health"',
        "immunity": 'immunity OR "immune system" OR "immune function" OR "immune response" OR "immune support" OR "immune health" OR "immune boosting" OR "immune enhancement" OR "immune modulation" OR "infection resistance"',
        "abdominal pain": 'abdominal pain OR "stomach pain" OR "belly pain" OR "abdominal discomfort" OR "gastrointestinal pain" OR "stomach ache" OR "gut pain" OR "visceral pain" OR "intestinal pain" OR "abdominal cramps" OR "stomach cramps" OR "epigastric pain" OR "gastrointestinal discomfort" OR "abdominal tenderness" OR "functional abdominal pain" OR "tummy ache" OR "tummy pain" OR "sore stomach" OR "sore abdomen" OR "abdomen pain" OR "bellyache" OR "gastralgia" OR "functional gut disorder" OR "indigestion pain" OR "dyspepsia" OR "gastric discomfort"',
        "digestion": 'digestion OR "digestive health" OR "digestive system" OR "digestive function" OR "digestive process" OR "gut health" OR "gastrointestinal health" OR "gastrointestinal function" OR "digestive enzymes" OR "gut function" OR "digestive issues" OR "digestive problems" OR "gut microbiome" OR "gut flora" OR "gut bacteria" OR "gut motility" OR "intestinal health" OR "digestive tract" OR "indigestion" OR "dyspepsia" OR "poor digestion" OR "slow digestion" OR "bloating" OR "gas" OR "flatulence" OR "upset stomach" OR "heartburn" OR "GERD" OR "acid reflux" OR "stomach upset" OR "leaky gut" OR "malabsorption" OR "digestive distress" OR "gastric emptying" OR "stomach acid" OR "bile"',
        "ibs": 'ibs OR "irritable bowel syndrome" OR "irritable colon" OR "spastic colon" OR "functional bowel disorder" OR "bowel irregularity" OR "intestinal disorder" OR "functional GI disorder" OR "altered bowel habits" OR "bowel dysfunction" OR "colonic spasm" OR "visceral hypersensitivity" OR "diarrhea-predominant IBS" OR "constipation-predominant IBS" OR "mixed IBS" OR "post-infectious IBS" OR "bowel issues" OR "bowel problems" OR "intestinal issues" OR "intestinal problems" OR "IBS-D" OR "IBS-C" OR "IBS-M" OR "sensitive bowel" OR "sensitive gut" OR "nervous bowel" OR "bowel pain" OR "intestinal pain" OR "functional bowel syndrome" OR "SIBO" OR "small intestinal bacterial overgrowth" OR "intestinal hypersensitivity" OR "bowel discomfort" OR "chronic bowel condition"'
    }
    
    # Check if any keys from our dictionary are in the query
    expanded_query_parts = []
    original_parts = query.lower().split()
    query_lower = query.lower()
    
    # First add the original query
    expanded_query_parts.append(f"({query})")
    
    # Check for matching terms and add their synonyms
    for term, synonyms in synonyms_dict.items():
        # Check if the key is in the query (including slight variations)
        term_variants = [term, term + "s", term + "es", term + "ing", term + "ed"]
        
        if any(variant in original_parts or any(variant in part for part in original_parts) for variant in term_variants):
            expanded_query_parts.append(f"({synonyms})")
            continue
            
        # Check if any of the synonyms are in the query
        # Extract individual synonyms from the value string
        synonym_list = [s.strip().lower().replace('"', '') for s in synonyms.split('OR')]
        
        # Check if any of these individual synonyms appear in the query
        if any(syn in query_lower or any(syn in part for part in original_parts) for syn in synonym_list):
            expanded_query_parts.append(f"({synonyms})")
    
    # If no matches found, just use the original query
    if len(expanded_query_parts) == 1:
        return query
        
    # Join all parts with AND to ensure both original query and expanded terms are included
    return " AND ".join(expanded_query_parts)

@cached_function
def fetch_abstracts(id_list):
    """Fetch abstracts for a list of PubMed IDs."""
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    return list(records)

def extract_abstract_info(records, original_query):
    """Extract relevant information from Medline records."""
    results = []
    
    # Process the original query to extract key terms
    query_terms = set(original_query.lower().split())
    # Get expanded query terms for context checking
    expanded_query_terms = get_expanded_terms_for_context(original_query)
    
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
            
            # Score the relevance of this abstract with context awareness
            relevance_score = score_relevance(abstract, title, supplements, original_query, expanded_query_terms)
            
            # Check if abstract is relevant enough (minimum threshold and must contain relevant context)
            # Higher threshold (5) to ensure real relevance
            if relevance_score >= 5 and is_context_relevant(abstract, title, expanded_query_terms):
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

def get_expanded_terms_for_context(query):
    """Get expanded terms for context relevance checking."""
    context_terms = set()
    query_lower = query.lower()
    
    # Health goal context mappings
    context_mappings = {
        "focus": ["focus", "attention", "concentration", "cognitive", "cognition", "brain", "mental", 
                 "memory", "alertness", "clarity", "performance", "executive function"],
        "energy": ["energy", "fatigue", "tired", "tiredness", "vitality", "stamina", "endurance", 
                  "lethargy", "exhaustion", "vigor"],
        "sleep": ["sleep", "insomnia", "rest", "circadian", "melatonin", "drowsiness", "wake", 
                 "rem", "deep sleep", "sleeping", "bedtime"],
        "stress": ["stress", "anxiety", "tension", "relaxation", "calm", "mood", "cortisol", 
                  "nervous", "anxious", "worry", "nervous system"],
        "inflammation": ["inflammation", "inflammatory", "anti-inflammatory", "cytokine", "immune", 
                        "pain", "chronic", "swelling"],
        "memory": ["memory", "recall", "cognitive", "brain", "cognition", "forgetfulness", 
                  "learning", "retention", "hippocampus", "dementia"],
        "abdominal pain": ["abdominal", "stomach", "belly", "pain", "ache", "discomfort", "cramps", 
                          "gastric", "intestinal", "visceral", "epigastric", "gi", "gastrointestinal",
                          "tummy", "abdomen", "bellyache", "sore", "gastralgia", "dyspepsia", 
                          "indigestion", "nausea", "tender", "spasm", "colic", "bloating"],
        "digestion": ["digestion", "digestive", "gut", "gastrointestinal", "intestinal", "bowel", 
                     "stomach", "microbiome", "enzymes", "flora", "bacteria", "motility", "gi", 
                     "gastric", "intestine", "colon", "indigestion", "dyspepsia", "bloating", "gas",
                     "flatulence", "heartburn", "gerd", "reflux", "malabsorption", "leaky gut",
                     "bile", "acid", "absorption", "metabolism", "transit", "empty", "digest", "stool"],
        "ibs": ["ibs", "irritable bowel", "bowel", "intestinal", "colon", "spastic", "diarrhea", 
               "constipation", "bloating", "gas", "abdominal", "functional disorder", "gi",
               "gut", "visceral", "hypersensitivity", "ibs-d", "ibs-c", "ibs-m", "altered habits",
               "sibo", "bacterial overgrowth", "sensitive bowel", "nervous bowel", "chronic bowel",
               "colonic", "flare", "flare-up", "trigger", "cramping", "urgency", "irregularity"]
    }
    
    # Add terms based on query content
    for key, terms in context_mappings.items():
        # Create common variations of the key for better matching
        key_variants = [key, key + "s", key + "es", key + "ing", key + "ed", 
                       "improve " + key, "enhance " + key, "better " + key, 
                       key + " issue", key + " problem", key + " concern",
                       key + " disorder", key + " condition", key + " health",
                       key + " support", key + " relief", key + " help"]
        
        # Check if any of the key variants are in the query
        if any(variant in query_lower for variant in key_variants):
            context_terms.update(terms)
    
    # Always add the original query terms
    context_terms.update(query_lower.split())
    
    return context_terms

def is_context_relevant(abstract, title, context_terms):
    """Check if the abstract's context is relevant to the user's query."""
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

@cached_function
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
    
    system_prompt = "You are a scientific research assistant specializing in nutritional supplements and health science."
    
    user_prompt = f"""Based on the following research abstracts about supplements for '{query}', provide a detailed analysis.
    
    Your analysis should address:
    
    1. Overview of the supplements mentioned and their potential benefits for the specific health goal
    2. Quality and strength of evidence for each supplement (clinical trials, animal studies, in vitro, etc.)
    3. Recommended dosages and safety considerations from the research
    4. Potential mechanisms of action explaining how these supplements work
    5. Note any conflicts between studies regarding effectiveness or recommended dosages
    
    Your report should be factual, evidence-based, and avoid exaggerating benefits or downplaying risks.
    If the evidence is inconclusive or limited, clearly state this fact.
    
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
    query = st.text_input("Health Goal", placeholder="e.g.'improve focus and mental clarity'")
    search_btn = st.button("Research Supplements", type="primary")

if search_btn:
    if not query:
        st.error("Please enter a health goal to research.")
    elif len(query.split()) < 4:
        st.warning("Be more specific. At least 4 words to describe your health goal for better results.")
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
                        
                        # Download button for the analysis
                        analysis_text = f"# Supplement Research Analysis for: {query}\n\n{analysis}\n\n## Sources\n\n"
                        for i, abstract in enumerate(abstracts, 1):
                            analysis_text += f"### {i}. {abstract['title']}\n"
                            analysis_text += f"*{abstract['authors']} - {abstract['journal']} ({abstract['year']})*\n\n"
                            analysis_text += f"{abstract['abstract']}\n\n"
                            analysis_text += f"PubMed ID: {abstract['pmid']}\n\n"
                        
                        st.download_button(
                            label="Download Analysis",
                            data=analysis_text,
                            file_name=f"supplement_research_{query.replace(' ', '_')}.md",
                            mime="text/markdown"
                        )

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