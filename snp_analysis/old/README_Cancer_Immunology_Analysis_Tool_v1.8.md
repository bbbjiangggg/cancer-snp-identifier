
# Cancer Immunology Analysis Tool v1.8

## Overview
This tool is designed to facilitate genomic analyses specifically tailored for cancer immunology research. It handles sequence data processing, mapping, and variant calling for given accession numbers of genomic data.

## Prerequisites
Before running this tool, ensure that you have Python 3.x installed along with the following packages:
- `pyfiglet`
- `termcolor`
- `subprocess`
- `os`
- `sys`
- `shutil`
- `socket`

Additionally, ensure that `bwa`, `bowtie2`, and `fastp` are installed and accessible from your system's PATH.

## Installation
1. Clone this repository or download the script directly to your local machine.
2. Install the required Python packages by running:
   ```
   pip install pyfiglet termcolor
   ```

## Configuration
1. Edit the script to include your specific paths to `bwa`, `bowtie2`, and `fastp` if they are not in the default locations (`/usr/local/bin` for `bwa` and `bowtie`, `/usr/bin` for `fastp`).
2. You need to configure your email settings within the script to enable the sending of notifications. Specify the SMTP details for your email provider.

## Usage
To run the script, navigate to the directory containing the script and type:
```
python cancer_immunology_analysis_v1.8.py
```
Follow the on-screen prompts to:
1. Enter your email for notifications.
2. Specify the job title for the analysis session.
3. Provide the path to the file containing accession numbers.
4. Indicate whether the reads are single-end or paired-end.
5. Enter specific chromosomes to analyze or type 'all' to analyze all chromosomes.
6. Specify how many accession numbers you wish to analyze.

## Input Data Format
The input file containing accession numbers should be a plain text file with one accession number per line.

## Output
The script will process each specified accession number and chromosome, generating a series of files including trimmed sequence files, SAM, BAM, and VCF files. Intermediate files will be deleted by default to save disk space.

## Troubleshooting
Ensure all paths are correctly specified and accessible. If errors occur, the script will output a red message indicating the failure reason. For interruptions, a notification will be sent to the specified email address.

## License
This script is distributed under the MIT license.

## Contact
For support or to report issues, please contact Robert Aguilar (robeagui@gmail.com).
