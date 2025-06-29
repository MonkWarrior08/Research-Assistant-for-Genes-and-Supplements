# Research Assistant for Genes and Supplements

This tool analyzes scientific literature from NCBI databases to provide insights on genes, SNPs, and supplements. It uses an AI model to summarize research papers and extract relevant information.

## Overview

This application is designed to help with researching genes, genetic variants (SNPs), and nutritional supplements. It retrieves scientific papers from NCBI, analyzes their content using an AI model, and presents the findings in a summarized format.

The goal is to make it easier to understand the implications of specific genetic markers and the scientific evidence behind various supplements.

## Features

- **Gene and SNP Research**: Search for scientific papers on specific genes, SNPs, and genotypes (e.g., "MTHFR C677T TT").
- **Supplement Research**: Look up studies related to nutritional supplements.
- **AI-Powered Analysis**: Get AI-generated summaries of research papers, focusing on key findings, health implications, and gene-supplement interactions.
- **NCBI Integration**: Directly queries PubMed and other NCBI databases.
- **Detailed Paper Information**: View abstracts, authors, journals, and publication dates for retrieved papers.
- **Two-Part Interface**: Separate pages for gene/SNP research and supplement research.

## Getting Started

### Prerequisites

- Python 3.8 or higher
- An OpenAI API key. An NCBI API key is optional but recommended.

### Installation

1.  **Clone the repository**
    ```bash
    git clone <repository-url>
    cd gene-research-assistant
    ```

2.  **Install dependencies**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Configure environment variables**
   
    Create a `.env` file in the root directory with the following:
    ```
    OPENAI_API_KEY=your_openai_api_key
    NCBI_EMAIL=your_email@example.com  # Required for NCBI E-utilities
    NCBI_API_KEY=your_ncbi_api_key     # Optional but recommended
    ```

### Running the Application

Launch the interactive web interface:
```bash
streamlit run Genes.py
```

## How to Use

### Gene Research

1.  **Enter a gene name, SNP ID, or combination**
    - Simple gene search: "MTHFR", "APOE", "COMT"
    - SNP search: "rs1801133", "rs429358"
    - Gene-SNP combination: "MTHFR rs1801133", "APOE rs429358"
    - With genotype: "MTHFR rs1801133 CT", "APOE rs429358 TT"

2.  **Customize your search (optional)**
    - Adjust the number of papers to retrieve
    - Filter by publication date
    - Include or exclude specific terms

3.  **Review research papers**
    - Expand paper entries to view detailed information
    - Sort results by relevance or date
    - Access direct links to original studies

4.  **Analyze genetic insights**
    - Review the AI-generated analysis of your genetic information
    - Learn about potential health implications
    - Discover lifestyle modifications relevant to your genetic profile
    - Understand how your genotype compares to population averages

### Supplement Research

1.  **Navigate to the Supplements page** using the sidebar
2.  **Enter a supplement name** to search for scientific research
3.  **Review the analysis** of supplement effectiveness and relevant studies

## Project Structure

- `Genes.py`: The main Streamlit application for gene research.
- `pages/supplements.py`: The page for supplement research.
- `ncbi_util.py`: Utility functions for interacting with NCBI databases.
- `openai_util.py`: Utility functions for OpenAI integration.
- `requirements.txt`: Project dependencies.
- `.env`: Environment variables (not tracked by git).

## Dependencies

This project uses several Python libraries, including:
- Biopython
- OpenAI
- Streamlit
- python-dotenv
- Requests
- BeautifulSoup4
- Pandas

See `requirements.txt` for the full list of dependencies.

## Privacy

This tool is designed for educational and informational purposes only. It should not be used to replace professional medical advice. Always consult with healthcare providers before making health decisions based on genetic information or supplement use.

The application does not store your search queries beyond your current session.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
