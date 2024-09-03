
# Genomic_Analysis_Tool_bam_v1.8

## Overview
This Python script is meticulously crafted to streamline the processing of genomic data, offering a robust solution for researchers in the field. It automates essential tasks such as downloading sequencing data from public databases, accurately mapping reads against reference genomes using BWA and Bowtie2, and efficiently cleaning up intermediate files to optimize storage. Additionally, the script features integration with SendGrid for real-time email notifications upon the successful completion of tasks or when interruptions occur. Notably, this script preserves the sorted BAM files, which are crucial for subsequent analyses like gene expression quantification, making it a comprehensive tool for genomic research.

## Prerequisites
Ensure Python 3.x is installed along with the following packages and utilities:
- `pyfiglet`
- `termcolor`
- `subprocess`
- `os`
- `sys`
- `shutil`
- `socket`
- `sendgrid`

You will also need:
- `bwa`
- `bowtie2`
- `fastp`

## Installation
1. Clone or download the script to your machine.
2. Install the required Python libraries:
   ```bash
   pip install pyfiglet termcolor sendgrid
   ```
3. Ensure the external utilities (`bwa`, `bowtie2`, `fastp`) are installed and properly configured in your PATH or specified directories.

## Configuration
Update the script with appropriate paths for:
- `bwa`
- `bowtie2`
- `fastp`

You must also configure the `SendGrid API` key in a `config.py` file like this:
```python
SENDGRID_API_KEY = 'your_sendgrid_api_key_here'
```

## Usage
Run the script from the command line:
```bash
python Genomic_Analysis_Tool_bam_v1.8
```
Follow the prompts to configure your analysis settings, including email for notifications, job title, path to the accession list file, read type (single-end or paired-end), and chromosomes to be analyzed.

## Input Data
Your accession list file should be a text file with one accession number per line.

## Outputs
The script outputs trimmed FASTQ files, SAM/BAM files, and optional VCF files. Intermediate files are deleted to conserve space.

## Email Notifications
Set up the tool to send notifications via email upon completion or interruption of tasks. Configure the sender email and ensure your SendGrid API key is valid.

## Troubleshooting
If errors occur, refer to the command line outputs colored in red. Make sure all file paths and configurations are correct. Email notifications will provide details on process interruptions.

## License
This script is freely distributed under the MIT license.

## Contact
For support or issues, contact Robert Aguilar (robeagui@gmail.com).
