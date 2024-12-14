import os
import pandas as pd
from termcolor import colored
import sys
import time
import threading
import re

skip_typing = False  # Global flag to indicate if typing should be skipped

def monitor_input():
    """Background thread that waits for the user to press Enter. Once pressed, sets skip_typing to True."""
    global skip_typing
    sys.stdin.readline()  # Wait until the user presses Enter
    skip_typing = True

# Start the background thread
input_thread = threading.Thread(target=monitor_input, daemon=True)
input_thread.start()

def live_type(text, color=None, delay=0.03, newline=True):
    """Print text to the screen as if typing it live. If user presses Enter, skip typing."""
    global skip_typing
    if color:
        text = colored(text, color)

    if skip_typing:
        # If skip is already activated, print instantly
        sys.stdout.write(text + ('\n' if newline else ''))
        sys.stdout.flush()
        return

    for i, char in enumerate(text):
        if skip_typing:
            # If skip was activated mid-typing, print the rest instantly
            remaining_text = text[i:]
            sys.stdout.write(remaining_text + ('\n' if newline else ''))
            sys.stdout.flush()
            return
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)

    if newline:
        sys.stdout.write('\n')
        sys.stdout.flush()

def parse_vcf(file_path):
    """Parses a VCF file and returns a DataFrame, handling rows with mismatched columns.
    If no data is found, returns an empty DataFrame. If the file doesn't exist, returns None."""
    if not os.path.exists(file_path):
        live_type(f"The file {file_path} does not exist.", "red")
        return None

    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if not line.startswith('##')]

    if not lines:
        live_type(f"No header or data lines found in {file_path}.", "yellow")
        # Return empty DataFrame indicating no data (not a parsing error)
        return pd.DataFrame()

    headers = lines[0].split("\t")
    data = []
    for line in lines[1:]:
        row = line.split("\t")
        if len(row) == len(headers):
            data.append(row)

    if data:
        df = pd.DataFrame(data, columns=headers)
        return df
    else:
        # Return empty DataFrame if no data rows were found after the header
        live_type(f"No data found in {file_path}.", "yellow")
        return pd.DataFrame(columns=headers)

def extract_cadd_score(info_field):
    """Extracts the CADD score from the INFO field."""
    cadd_score = None
    if 'CADD' in info_field:
        try:
            cadd_info = [item for item in info_field.split(';') if item.startswith('CADD=')]
            if cadd_info:
                cadd_value = cadd_info[0].split('=')[1].split(':')[-1]
                cadd_score = float(cadd_value)
        except (IndexError, ValueError):
            cadd_score = None
    return cadd_score

def extract_cosmic_info(info_field):
    """Extracts COSMIC information from the INFO field."""
    cosmic_info = None
    if 'COSMIC' in info_field:
        try:
            cosmic_data = [item for item in info_field.split(';') if item.startswith('COSMIC=')]
            if cosmic_data:
                cosmic_info = cosmic_data[0].split('=')[1]
        except IndexError:
            cosmic_info = None
    return cosmic_info

def extract_important_snps():
    vcf_files = [
        'cpg.vcf', 'ensembl.vcf', 'cadd.vcf', 'clinvar.vcf', 
        'cosmic.vcf', 'gnomad.vcf', 'gwas.vcf', 'phast.vcf', 'sift.vcf'
    ]

    important_snps = {}

    live_type("Press Enter at any time to skip the live typing effect.", "magenta")
    live_type("Extracting important SNPs from the following VCF files:", "cyan")
    for f in vcf_files:
        live_type(f" - {f}", "cyan", delay=0.01)

    for file_name in vcf_files:
        if file_name.endswith('.vcf'):
            if os.path.exists(file_name):
                live_type(f"Parsing {file_name}...", "blue")
                df = parse_vcf(file_name)
                if df is None:
                    # df is None could mean file not found or a more severe parsing issue
                    # already handled above by message printing.
                    live_type(f"Skipping {file_name} due to parsing issues (file not found or other error).", "red")
                else:
                    # df is a DataFrame (possibly empty)
                    if df.empty:
                        # It's not a parsing error, just no data
                        live_type(f"No data in {file_name}, skipping processing.", "yellow")
                    else:
                        # Proceed with filtering if not empty
                        if 'INFO' in df.columns:
                            live_type(f"Filtering important SNPs in {file_name}...", "yellow")
                            for index, row in df.iterrows():
                                info_field = row['INFO']
                                cadd_score = extract_cadd_score(info_field)
                                cosmic_info = extract_cosmic_info(info_field)
                                snp_id = row['ID']
                                position = row['POS']
                                chromosome = row['#CHROM']

                                # Unique key for the SNP
                                snp_key = (chromosome, position, snp_id)

                                # Criteria check
                                if (
                                    'COSMIC' in info_field or 'ClinVar' in info_field or
                                    (cadd_score is not None and cadd_score > 20) or
                                    ('SIFT' in info_field and 'deleterious' in info_field) or
                                    'GWAS' in info_field
                                ):
                                    if snp_key not in important_snps:
                                        important_snps[snp_key] = {
                                            'Chromosome': chromosome,
                                            'Position': position,
                                            'ID': snp_id,
                                            'Reference': row['REF'],
                                            'Alternate': row['ALT'],
                                            'CADD Score': cadd_score if cadd_score is not None else 'N/A',
                                            'COSMIC Info': cosmic_info if cosmic_info is not None else 'N/A',
                                            'INFO': info_field,
                                            'Source File': file_name
                                        }
                                    else:
                                        # Update COSMIC info if found again
                                        if cosmic_info is not None:
                                            important_snps[snp_key]['COSMIC Info'] = cosmic_info
                        else:
                            live_type(f"File {file_name} does not contain 'INFO' column. Skipping.", "red")
            else:
                live_type(f"File {file_name} does not exist. Skipping.", "red")
        else:
            live_type(f"{file_name} is not a VCF file. Skipping.", "red")

    # Save results
    if important_snps:
        important_snps_df = pd.DataFrame.from_dict(important_snps, orient='index')
        important_snps_df.to_csv('summary_of_important_snps.csv', index=False)
        live_type("Important SNPs saved to 'summary_of_important_snps.csv'", "green")
    else:
        no_data_df = pd.DataFrame([["No important SNPs found based on the criteria."]], columns=["Message"])
        no_data_df.to_csv('summary_of_important_snps.csv', index=False)
        live_type("No important SNPs found based on the criteria. CSV with message generated.", "yellow")

if __name__ == "__main__":
    extract_important_snps()
