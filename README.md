# Gene Research Assistant

This application helps researchers analyze scientific literature related to specific genes. It:

1. Takes a gene name/ID as input
2. Queries the NCBI databases via their E-utilities API to find relevant research papers
3. Uses OpenAI to generate insights and analysis based on the retrieved papers

## Setup

1. Clone this repository
2. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```
3. Create a `.env` file with your API keys:
   ```
   OPENAI_API_KEY=your_openai_api_key
   NCBI_API_KEY=your_ncbi_api_key  # Optional but recommended for higher rate limits
   ```

## Usage

Run the Streamlit app:
```
streamlit run app.py
```

## Features

- Search for research papers by gene name
- Filter results by publication date, relevance, etc.
- Generate AI-powered summaries and insights
- View key findings across multiple papers

## Technical Implementation

- Uses Biopython's Entrez API to access NCBI databases
- Integrates with OpenAI's API for analysis
- Built with Streamlit for a simple, interactive interface 