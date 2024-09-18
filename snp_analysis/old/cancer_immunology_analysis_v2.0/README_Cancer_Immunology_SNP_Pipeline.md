
# Cancer Immunology v2.0 - SNP Analysis Pipeline

## Description
This script automates the processing and analysis of SNP data for cancer immunology research. It performs the following tasks:
- Checks and installs required Python packages.
- Downloads and trims sequence data.
- Maps the sequence data to reference genomes using BWA and Bowtie2.
- Performs variant calling to generate VCF files.
- Logs all actions with color-coded messages for easy identification of warnings, errors, and success messages.

## Features
- **Package Management**: Automatically checks and installs missing Python packages.
- **File Management**: Ensures directories exist, handles file deletions, and clears directories as needed.
- **Chromosome Handling**: Analyzes specified chromosomes, with the option to process all chromosomes.
- **Logging**: Provides detailed, color-coded logs to track the process and identify any issues.
- **Multi-threading**: Utilizes multiple CPU threads to speed up the analysis.

## Prerequisites
- **Python 3.x**
- **BWA**: Ensure BWA is installed and accessible via the specified path.
- **Bowtie2**: Ensure Bowtie2 is installed and accessible via the specified path.
- **Samtools**: Required for handling SAM/BAM files.
- **BCFtools**: Required for variant calling and generating VCF files.
- **fastp**: Required for trimming sequence data.

## Installation
1. Clone the repository or download the script file.
2. Ensure all dependencies (BWA, Bowtie2, Samtools, BCFtools, and fastp) are installed on your system.
3. Run the script using Python:
   ```bash
   python3 script_name.py
   ```

## Usage

1. **BWA and Bowtie Paths**: Ensure that the base paths for BWA and Bowtie are correctly set or entered when prompted.
2. **Accession List**: Provide the path to a file containing accession numbers (one per line) when prompted.
3. **Read Type**: Indicate whether the reads are single-end (1) or paired-end (2).
4. **Threads**: Specify the number of CPU threads to use for the analysis.
5. **Chromosomes**: Specify the chromosomes to analyze (e.g., `1,2,3`), or type `all` to process all chromosomes.
6. **Run the Analysis**: The script will process each accession number, download the sequences, trim them, and perform mapping and variant calling.

### Example
```bash
python3 script_name.py
```

When prompted:
- Enter the accession list file path.
- Choose between single-end or paired-end reads.
- Specify the number of threads (e.g., 4).
- Enter the chromosomes to analyze (e.g., `1,2,X,Y`) or `all`.

## Logging
- **Info**: General information, displayed in magenta.
- **Warning**: Alerts for issues that won't stop the script, displayed in yellow.
- **Error**: Critical errors that cause the script to stop, displayed in red.
- **Success**: Successful operations, displayed in green.

## Author
Dr. Aguilar

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
