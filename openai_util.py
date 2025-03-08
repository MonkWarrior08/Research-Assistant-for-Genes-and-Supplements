import os
from typing import Dict, List, Any, Optional
import openai
from dotenv import load_dotenv

load_dotenv()

openai.api_key - os.getenv("OPENAI_API_KEY")
if not openai.api_key:
    print("openai api key is not set")