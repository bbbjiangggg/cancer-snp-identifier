import os
from datetime import datetime
import pandas as pd
from io import StringIO

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Function to convert VCF to CSV using pandas
def vcf_to_csv(vcf_file, csv_file):
    with open(vcf_file, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(
        StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )
    df = df[['#CHROM', 'POS', 'REF', 'ALT']]
    df.columns = ['Chrom', 'Pos', 'Ref', 'Mut']
    df.to_csv(csv_file, index=False)

# Set the parent directory path
parent_directory_path = os.getcwd()
print(f"{MAGENTA}\nParent directory path:{RESET} {parent_directory_path} ")

# Prompt the user for the name of the VCF file they want to convert to CSV
vcf_file_name = input(f'{MAGENTA}\nEnter the name of the VCF file you want to convert to CSV:{RESET} ')
vcf_file_path = os.path.join(parent_directory_path, vcf_file_name)

# Convert the VCF file to CSV
csv_file_name = f'{vcf_file_name}.csv'
csv_file_path = os.path.join(parent_directory_path, csv_file_name)
print(f'{MAGENTA}\nConverting {vcf_file_name} to CSV...{RESET}')
vcf_to_csv(vcf_file_path, csv_file_path)

if os.path.isfile(csv_file_path):
    print(f'{MAGENTA}Conversion complete. CSV file saved at:{RESET} {csv_file_path}')

