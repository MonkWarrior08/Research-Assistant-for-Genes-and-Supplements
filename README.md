# 🧬 Gene Research Assistant

![Version](https://img.shields.io/badge/version-1.0.0-blue)
![Python](https://img.shields.io/badge/Python-3.8+-green.svg)
![License](https://img.shields.io/badge/license-MIT-yellow.svg)

> **Unlock your genetic insights** - A powerful tool for analyzing scientific literature related to specific genes and SNPs to provide personalized health insights.

<img width="1760" alt="Screenshot 2025-03-08 at 9 07 24 PM" src="https://github.com/user-attachments/assets/5a41d589-5c45-4a0a-bf9a-efc6116f9846" /><img width="1748" alt="Screenshot 2025-03-08 at 9 07 55 PM" src="https://github.com/user-attachments/assets/7391e5e2-56af-4877-8192-a417b668fd90" />



## 📋 Overview

Gene Research Assistant is an AI-powered application that helps individuals understand scientific literature related to specific genes and genetic variants (SNPs). The tool bridges the gap between complex scientific research and personal genetic insights by:

- Retrieving relevant scientific papers from NCBI databases
- Processing and analyzing paper content with AI
- Providing personalized, accessible interpretations of genetic information
- Translating complex scientific findings into practical health insights

## ✨ Features

### 🔍 Research Paper Analysis
- **Gene-specific searches**: Find papers related to a particular gene or SNP
- **Comprehensive paper details**: View titles, authors, journals, publication dates, and full abstracts
- **Efficient filtering**: Organize research findings by relevance

### 🧠 AI-Powered Genetic Insights
- **Personal genetic insights**: Understand what your genes mean for your health and traits
- **Practical applications**: Get actionable lifestyle recommendations based on genetic information
- **SNP profile analysis**: Learn specific implications of your genetic variants
- **Plain language explanations**: Complex scientific concepts explained in accessible terms

### 🛠️ Technical Capabilities
- **NCBI Integration**: Direct access to PubMed and other scientific databases
- **Real-time analysis**: Generate insights on-demand from the latest research
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
streamlit run app.py
```

## 📖 How to Use

1. **Enter a gene name or SNP ID**
   - Example formats: "BRCA1", "rs429358", "APOE rs429358"
   - You can also include genotype information: "APOE rs429358 TT" 

2. **Review research papers**
   - Expand paper entries to view detailed information
   - Read abstracts to understand key findings

3. **Analyze genetic insights**
   - View the AI-generated analysis for personalized information
   - Learn about practical applications for your health

## 🧩 Project Structure

- `app.py` - Main Streamlit application and user interface
- `ncbi_util.py` - Functions for interacting with NCBI databases
- `openai_util.py` - OpenAI integration for AI-powered analysis
- `requirements.txt` - Project dependencies
- `.env` - Environment variables and API keys (not tracked by git)

## 📚 Dependencies

- **Biopython**: Interface with NCBI's Entrez databases
- **OpenAI**: AI-powered analysis of scientific literature
- **Streamlit**: Interactive web interface
- **Python-dotenv**: Environment variable management
- **Requests**: HTTP requests handling

## 🔒 Privacy & Ethics

This tool is designed for educational and informational purposes only. It should not be used to replace professional medical advice. Always consult with healthcare providers before making health decisions based on genetic information.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

<p align="center">Made with ❤️ for genetic discovery and personal health insights</p> 
