import streamlit as st
import time
import os
from dotenv import load_dotenv
import re
import openai

from ncbi_util import search_gene_paper, fetch_paper_details, get_gene_info
from openai_util import analyze_papers, analyze_single_abstract

load_dotenv()

st.set_page_config(
    page_title="Research Assistant",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Add a sidebar menu for navigation (this is optional since Streamlit will automatically
# add navigation for pages in the pages/ directory)
with st.sidebar:    
    st.header("Search Parameters")
    st.markdown("Enter gene information:")
    gene_symbol = st.text_input("Gene Symbol", value="", placeholder="e.g. OPN4")
    rs_id = st.text_input("SNP Identifier", value="", placeholder="e.g. rs1079610")
    genotype = st.text_input("Genotype", value="", placeholder="e.g. TC")
    search_btn = st.button("Analyze", type='primary')
    
    st.markdown("---")
    st.header("About Gene Research")
    st.markdown("""
    This assistant uses:
    - **NCBI/PubMed API** for scientific gene research and paper retrieval
    - **OpenAI GPT-4o** for intelligent analysis of research findings
    """)

st.title("Gene Research Assistant")
st.markdown("""
This AI assistant searches the National Center for Biotechnology Information (NCBI) 
for scientific research on genes and provides a detailed analysis.

Enter a gene name or SNP identifier in the sidebar to begin.
""")

if search_btn:
    if not rs_id or not genotype:
        st.error("Please fill at least SNP Identifier and Genotype")
    elif not re.search(r'rs\d{3,}', rs_id.lower()):
        st.error("SNP Identifier must be in the format 'rs' followed by at least 3 digits (e.g., rs1079610)")
    elif not re.match(r'^[ACGTacgt]{2}$', genotype):
        st.error("Genotype must be exactly two nucleotides (A, C, G, or T)")
    else:
        # Convert inputs to consistent format (gene_symbol as is, rs_id lowercase, genotype uppercase)
        gene_symbol = gene_symbol.upper()
        rs_id = rs_id.lower()
        genotype = genotype.upper()
        # Combine inputs for processing
        gene_name = f"{gene_symbol} {rs_id} {genotype}" if gene_symbol else f"{rs_id} {genotype}"
        
        tabs = st.tabs(["Research Papers", "Analysis"])
        
        with st.spinner(f"Researching for {gene_name}"):
            combined_info = get_gene_info(gene_name)
            max_results = 20
            custom_date_range = None
            
            # Create a list to collect all individual paper analyses
            all_analyses = []

            with tabs[0]:
                st.header(f"Research papers for {gene_name}")

                pmids = search_gene_paper(gene_name, max_results, custom_date_range)

                if not pmids:
                    st.warning(f"No papers found for {gene_name}. Try different parameters")
                else:
                    st.success(f"Found {len(pmids)} papers.")
                    
                    papers = fetch_paper_details(pmids)

                    for i, paper in enumerate(papers):
                        with st.expander(f"{i+1}. {paper['title']}"):
                            st.markdown(f"**Authors:** {paper['author']}")
                            st.markdown(f"**Journal:** {paper['journal']}")
                            st.markdown(f"**Publication Date:** {paper['publication_date']}")
                            if paper['doi']:
                                st.markdown(f"**DOI:** {paper['doi']}")
                            st.markdown(f"**PMID:** {paper['pmid']}")
                            st.markdown(f"**URL:** [PubMed Link]({paper['url']})")
                            
                            if paper['abstract']:
                                # Display full abstract first
                                st.markdown("### Full Abstract")
                                st.markdown(paper['abstract'])
                                
                                # Generate and display AI analysis for this specific abstract
                                with st.spinner("Analyzing this paper..."):
                                    analysis_result = analyze_single_abstract(paper['abstract'], gene_name, rs_id, genotype)
                                    
                                    # Store the analysis for the summary tab if it contains relevant information
                                    if "No relevant information found" not in analysis_result:
                                        paper_info = {
                                            "title": paper['title'],
                                            "author": paper['author'],
                                            "journal": paper['journal'],
                                            "publication_date": paper['publication_date'],
                                            "analysis": analysis_result
                                        }
                                        all_analyses.append(paper_info)
                                    
                                    # Show AI analysis
                                    st.markdown("### AI Analysis")
                                    st.markdown("*The following analysis was automatically generated for this research paper:*")
                                    
                                    if "No relevant information found" in analysis_result:
                                        st.info(analysis_result)
                                    else:
                                        st.success(analysis_result)
                            else:
                                st.info("No abstract available for this paper.")
            
            # Analysis tab - Summarize all individual analyses
            with tabs[1]:
                st.header(f"Summary Analysis for {gene_name}")
                
                if not all_analyses:
                    st.warning("No relevant information found in any of the papers to create a summary analysis.")
                else:
                    with st.spinner("Generating comprehensive summary analysis..."):
                        # Format the analyses for the summary
                        analyses_text = ""
                        for i, paper_info in enumerate(all_analyses):
                            analyses_text += f"## PAPER {i+1}: {paper_info['title']}\n"
                            analyses_text += f"*{paper_info['author']} - {paper_info['journal']} ({paper_info['publication_date']})*\n\n"
                            analyses_text += f"{paper_info['analysis']}\n\n"
                        
                        # Generate a comprehensive summary of all analyses
                        summary_prompt = f"""
                        Based on the following individual paper analyses about {gene_name}, create a comprehensive summary.
                        These are analyses of papers about the gene {gene_symbol.upper() if gene_symbol else ''}, 
                        SNP {rs_id}, and genotype {genotype}.
                        
                        Focus on synthesizing the key findings and creating a cohesive summary that highlights:
                        1. The most important findings about the SNP and genotype across all papers
                        2. Any consensus or contradictions between the papers
                        3. The most significant health implications and practical considerations
                        
                        Format your response with clear sections and bullet points for readability.
                        
                        INDIVIDUAL PAPER ANALYSES:
                        {analyses_text}
                        """
                        
                        response = openai.chat.completions.create(
                            model="gpt-4o",
                            messages=[
                                {"role": "system", "content": "You are a scientific research assistant specializing in genetics. You synthesize multiple research findings into coherent summaries."},
                                {"role": "user", "content": summary_prompt}
                            ],
                            temperature=0.3,
                            max_tokens=2000
                        )
                        
                        summary_analysis = response.choices[0].message.content.strip()
                        st.markdown(summary_analysis)
                        
                        # Add download button for the summary analysis
                        st.download_button(
                            label="Download summary analysis",
                            data=summary_analysis,
                            file_name=f"{gene_symbol or rs_id}_{genotype}_summary_analysis.txt",
                            mime="text/plain"
                        )

st.markdown("---")

if not os.getenv("OPENAI_API_KEY"):
    st.sidebar.warning("Openai api key not found. Make sure it's added in the .env file")