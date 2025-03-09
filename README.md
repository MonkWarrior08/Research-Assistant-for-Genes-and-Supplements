# 🧬 Research Assistant for Genes and Supplements

![Version](https://img.shields.io/badge/version-1.1.0-blue)
![Python](https://img.shields.io/badge/Python-3.8+-green.svg)
![License](https://img.shields.io/badge/license-MIT-yellow.svg)

> **Unlock your health insights** - A powerful tool for analyzing scientific literature related to specific genes, SNPs, and supplements to provide personalized health recommendations.


## 📋 Overview

The Research Assistant for Genes and Supplements is an AI-powered application that helps individuals understand scientific literature related to specific genes, genetic variants (SNPs), and nutritional supplements. The tool bridges the gap between complex scientific research and personal health insights by:

- Retrieving relevant scientific papers from NCBI databases
- Processing and analyzing paper content with AI
- Providing personalized, accessible interpretations of genetic information
- Translating complex scientific findings into practical health insights
- Offering supplement research based on scientific literature

## ✨ Features

### 🔍 Research Paper Analysis
- **Gene-specific searches**: Find papers related to a particular gene or SNP
- **Supplement research**: Analyze scientific literature on nutritional supplements
- **Comprehensive paper details**: View titles, authors, journals, publication dates, and full abstracts
- **Efficient filtering**: Organize research findings by relevance

### 🧠 AI-Powered Insights
- **Personal genetic insights**: Understand what your genes mean for your health and traits
- **Practical applications**: Get actionable lifestyle recommendations based on genetic information
- **SNP profile analysis**: Learn specific implications of your genetic variants
- **Supplement effectiveness**: Research evidence-based information on supplements
- **Plain language explanations**: Complex scientific concepts explained in accessible terms

### 🛠️ Technical Capabilities
- **NCBI Integration**: Direct access to PubMed and other scientific databases
- **Real-time analysis**: Generate insights on-demand from the latest research
- **Multi-page interface**: Separate sections for gene and supplement research
- **User-friendly interface**: Simple, intuitive design for both researchers and individuals

## 🚀 Getting Started

### Prerequisites
- Python 3.8 or higher
- API keys for OpenAI and NCBI (optional for NCBI, but recommended)

### Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd gene-research-assistant
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Configure environment variables**
   
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

## 📖 How to Use

### Gene Research
1. **Enter a gene name or SNP ID**
   - Example formats: "BRCA1", "rs429358", "APOE rs429358"
   - You can also include genotype information: "APOE rs429358 TT" 

2. **Review research papers**
   - Expand paper entries to view detailed information
   - Read abstracts to understand key findings

3. **Analyze genetic insights**
   - View the AI-generated analysis for personalized information
   - Learn about practical applications for your health

### Supplement Research
1. **Navigate to the Supplements page** using the sidebar
2. **Enter a supplement name** to search for scientific research
3. **Review the analysis** of supplement effectiveness and relevant studies

## 💻 Project Structure

- `Genes.py` - Main Streamlit application for gene research
- `pages/supplements.py` - Supplement research functionality
- `ncbi_util.py` - Functions for interacting with NCBI databases
- `openai_util.py` - OpenAI integration for AI-powered analysis
- `cache_utils.py` - Caching system for improved performance
- `requirements.txt` - Project dependencies
- `.env` - Environment variables and API keys (not tracked by git)

## 📚 Dependencies

- **Biopython**: Interface with NCBI's Entrez databases
- **OpenAI**: AI-powered analysis of scientific literature
- **Streamlit**: Interactive web interface
- **Python-dotenv**: Environment variable management
- **Requests**: HTTP requests handling
- **BeautifulSoup**: HTML parsing for supplementary data (used in supplements research)

## 🔒 Privacy & Ethics

This tool is designed for educational and informational purposes only. It should not be used to replace professional medical advice. Always consult with healthcare providers before making health decisions based on genetic information or supplement use.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

<p align="center">Made with ❤️ for scientific discovery and personal health insights</p> 
