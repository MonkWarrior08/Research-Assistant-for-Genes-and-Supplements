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

