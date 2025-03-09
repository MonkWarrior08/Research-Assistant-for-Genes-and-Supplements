import os
import streamlit as st
import openai
from Bio import Entrez, Medline
import requests
from bs4 import BeautifulSoup
import pandas as pd
from dotenv import load_dotenv
import time
import re
import json
from cache_utils import cached_function

load_dotenv()

openai.api_key = os.getenv("OPENAI_API_KEY")

Entrez.email = os.getenv("ENTREZ_EMAIL")

st.set_page_config(
    page_title="Supplement Research Assistant",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

def search_pubmed(query, max_result=20):
    expanded_term = expand_query(query)
    search_term = f"({expanded_term}) AND (supplement OR supplementation OR nutraceutical OR herbal OR dietary OR natural compound OR natural product OR vitamin OR mineral OR nootropic OR cognitive enhancer)"
    
    handle = Entrez.esearch(
        db="pubmed",
        term=search_term,
        retmax=max_result,
        sort="relevance",
        field="title/abstract"
    )
    results = Entrez.read(handle)
    handle.close()
    return results["IdList"]

def expand_query(query):
    synonyms_dict = {
        "focus": 'focus OR attention OR concentration OR cognition OR cognitive OR "mental clarity" OR alertness OR "mental performance" OR "brain function" OR "cognitive performance" OR "cognitive enhancement" OR "mental focus" OR "attention span" OR "cognitive ability"',
        "energy": 'energy OR fatigue OR tiredness OR "mental energy" OR vitality OR stamina OR endurance OR "physical performance" OR "mental fatigue" OR "physical fatigue" OR "energy boost" OR "energy level"',
        "sleep": 'sleep OR insomnia OR "sleep quality" OR "sleep duration" OR "sleep disorder" OR "sleep disturbance" OR "sleep onset" OR "sleep maintenance" OR "deep sleep" OR "REM sleep" OR "sleep efficiency" OR "sleep architecture"',
        "stress": 'stress OR anxiety OR tension OR "stress relief" OR "stress reduction" OR relaxation OR "cortisol level" OR "stress hormone" OR "stress management" OR "anxiety reduction" OR calm OR "mental stress" OR "chronic stress"',
        "inflammation": 'inflammation OR "inflammatory marker" OR "anti-inflammatory" OR "chronic inflammation" OR cytokine OR "inflammatory response" OR "systemic inflammation" OR "inflammatory cytokine" OR "inflammation reduction"',
        "memory": 'memory OR "memory enhancement" OR "memory formation" OR "memory recall" OR "memory retention" OR "short-term memory" OR "long-term memory" OR "working memory" OR "memory function" OR "memory improvement" OR "memory consolidation"',
        "mood": 'mood OR depression OR "mood enhancement" OR "mood improvement" OR "mood regulation" OR "emotional wellbeing" OR "positive mood" OR "mood stabilization" OR "mood disorder" OR "mood support" OR happiness OR "emotional health"',
        "immunity": 'immunity OR "immune system" OR "immune function" OR "immune response" OR "immune support" OR "immune health" OR "immune boosting" OR "immune enhancement" OR "immune modulation" OR "infection resistance"'
    }
    expanded_query_parts = []
    original_parts = query.lower().split()

    expanded_query_parts.append(f"({query})")

    for term, synonyms in synonyms_dict.items():
        if term in original_parts or any(term in part for part in original_parts):
            expanded_query_parts.append(f"({synonyms})")
    
    if len(expanded_query_parts) == 1:
        return query
    
    return " OR ".join(expanded_query_parts)

@cached_function
def fetch_abstract(id_list):
    ids = ",".join(id_list)
    handle = Entrez.efetch(
        db="pubmed",
        id=ids,
        rettype="medline",
        retmode="text"
    )
    records = Medline.parse(handle)
    return list(records)