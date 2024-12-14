import os
import re
import pandas as pd
import glob
import shutil
import requests
from termcolor import colored
import sys
import time
import threading

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
    """Print text to the screen as if typing it live. If user presses Enter at any time, skip typing."""
    global skip_typing
    if color:
        text = colored(text, color)

    if skip_typing:
        # If skip is already activated, print instantly
        sys.stdout.write(text + ('\n' if newline else ''))
        sys.stdout.flush()
        return

    for char in text:
        if skip_typing:
            # If skip was activated mid-typing, print the rest instantly
            sys.stdout.write(text[text.index(char):] + ('\n' if newline else ''))
            sys.stdout.flush()
            return
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)

    if newline:
        sys.stdout.write('\n')
        sys.stdout.flush()

def parse_position(position_str):
    """Convert a chromosome position string with commas to an integer."""
    return int(re.sub(r'[^\d]', '', position_str))

def parse_position_range(range_str):
    """Parse a range of chromosome positions formatted as 'start-end'."""
    try:
        start_str, end_str = range_str.split('-')
        start_position = parse_position(start_str)
        end_position = parse_position(end_str)
        return start_position, end_position
    except ValueError:
        live_type("Invalid range format. Please enter the chromosome position range as 'start-end'.", "red")
        return None, None

def query_gene_positions(gene_name):
    """Query Ensembl for the positions of a given gene and return the chromosome, start, and end positions."""
    url = f"https://rest.ensembl.org/lookup/symbol/human/{gene_name}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        live_type("Gene: " + data['display_name'], "cyan")
        live_type("Chromosome: " + str(data['seq_region_name']), "cyan")
        live_type("Start Position: " + str(data['start']), "cyan")
        live_type("End Position: " + str(data['end']), "cyan")
        return data['seq_region_name'], data['start'], data['end']
    else:
        live_type(f"Error fetching data for gene '{gene_name}': {response.status_code}", "red")
        if response.text:
            live_type(f"Response: {response.text}", "red")
        return None, None, None

def find_vcf_file(chromosome):
    """Find the VCF file in the directory corresponding to the given chromosome."""
    chromosome_dir = f"ch{chromosome}"

    # Check if the chromosome directory exists
    if not os.path.isdir(chromosome_dir):
        live_type(f"Directory '{chromosome_dir}' not found. Please ensure the correct directory exists.", "red")
        return None

    os.chdir(chromosome_dir)
    live_type(f"Entered directory: {chromosome_dir}", "blue")

    # Pattern: Accept letters and/or numbers for the cancer type field
    pattern = re.compile(rf"ch{chromosome}_[0-9A-Za-z]+_comb\.vcf", re.IGNORECASE)
    vcf_files = [f for f in os.listdir('.') if pattern.match(f)]

    if vcf_files:
        live_type("Detected VCF files:", "yellow")
        for idx, file in enumerate(vcf_files, 1):
            live_type(f"{idx}. {file}", "yellow", delay=0.01)

        if len(vcf_files) == 1:
            selected_file = vcf_files[0]
        else:
            choice = input(colored("Please enter the number of the correct VCF file, or '0' if none: ", "magenta"))
            if choice == '0':
                live_type("No VCF file selected. Exiting.", "red")
                return None
            try:
                selected_file = vcf_files[int(choice) - 1]
            except (ValueError, IndexError):
                live_type("Invalid choice. Exiting.", "red")
                return None

        live_type(f"Found VCF file: {selected_file}", "green")
        return selected_file
    else:
        live_type("No VCF files found matching the pattern.", "red")
        return None

def filter_positions():
    # Display initial instructions to the user with live typing
    live_type("Press Enter at any time to skip the live typing effect.", "magenta")
    live_type("Welcome to the SNP Filtering and Processing Tool!", "magenta")
    live_type("This program will:", "cyan")
    live_type("- Take a gene name and use Ensembl to find its chromosome coordinates.", "cyan")
    live_type("- Locate a VCF file in a directory corresponding to that chromosome (e.g., ch22).", "cyan")
    live_type("- The VCF file must follow the naming convention: ch<chromosome>_<cancer_type>_comb.vcf", "cyan")
    live_type("   where <cancer_type> can be letters and/or numbers.", "cyan")
    live_type("Please ensure that:", "yellow")
    live_type("1) A directory named 'ch<chromosome_number>' already exists.", "yellow")
    live_type("2) Inside it, there is a combined SNP VCF file matching the pattern (e.g., ch22_ovca2_comb.vcf).", "yellow")
    live_type("After processing, this script will filter variants based on gene coordinates, retrieve SNP IDs, and output a processed CSV and a list of SNPs.", "cyan")
    live_type("", newline=True)  # Just a blank line for spacing

    gene_name = input(colored("Please enter the gene name: ", "magenta"))
    chromosome, start_position, end_position = query_gene_positions(gene_name)

    if chromosome is None or start_position is None or end_position is None:
        # If gene information not found
        live_type("No valid gene information retrieved. Please check the gene name and try again.", "red")
        return

    # Automatically navigate to the chromosome directory and find VCF file
    input_file = find_vcf_file(chromosome)
    if not input_file:
        return  # Exit if no file is found or selected

    # Create a directory for the gene if it doesn't exist
    gene_dir = os.path.join(os.getcwd(), gene_name)
    os.makedirs(gene_dir, exist_ok=True)
    live_type(f"Results will be saved in the directory: {gene_dir}", "blue")

    # Determine if the file is a CSV or VCF based on the file extension
    file_extension = input_file.split('.')[-1].lower()

    if file_extension == 'csv':
        try:
            df = pd.read_csv(input_file)
        except FileNotFoundError:
            live_type("The file was not found. Please check the path and try again.", "red")
            return
        except pd.errors.EmptyDataError:
            live_type("The file is empty. Please provide a valid CSV file.", "red")
            return
        except Exception as e:
            live_type(f"An error occurred: {e}", "red")
            return

        required_columns = {'Chromosome', 'Position'}
        if not required_columns.issubset(df.columns):
            live_type("The CSV file does not contain the required headers ('Chromosome', 'Position'). Please check the file format.", "red")
            return

        df['Position'] = df['Position'].apply(lambda x: parse_position(str(x)))
        filtered_df = df[(df['Position'] >= start_position) & (df['Position'] <= end_position)].copy()

        output_file = os.path.join(gene_dir, "filtered_positions.csv")
        filtered_df.to_csv(output_file, index=False)
        live_type(f"Filtered data saved to {output_file}", "green")

    elif file_extension == 'vcf':
        try:
            df = pd.read_csv(input_file, delim_whitespace=True, header=None, usecols=[0, 1, 4])
        except FileNotFoundError:
            live_type("The file was not found. Please check the path and try again.", "red")
            return
        except pd.errors.EmptyDataError:
            live_type("The file is empty. Please provide a valid VCF file.", "red")
            return
        except Exception as e:
            live_type(f"An error occurred: {e}", "red")
            return

        df.columns = ['Chromosome', 'Position', 'Genotype']
        df['Position'] = df['Position'].apply(lambda x: parse_position(str(x)))
        filtered_df = df[(df['Position'] >= start_position) & (df['Position'] <= end_position)].copy()

        def calculate_ratio(genotype):
            count_ones = genotype.count('1')
            count_zeros = genotype.count('0')
            total = count_ones + count_zeros
            return count_ones / total if total > 0 else 0

        filtered_df.loc[:, 'Ratio_1s'] = filtered_df['Genotype'].apply(calculate_ratio)
        filtered_df = filtered_df.drop(columns=['Genotype'])
        filtered_df = filtered_df.sort_values(by='Ratio_1s', ascending=False)

        output_file = os.path.join(gene_dir, "filtered_positions_vcf.csv")
        filtered_df.to_csv(output_file, index=False)
        live_type(f"Filtered data saved to {output_file}", "green")

    else:
        live_type("Unsupported file type. Please provide a CSV or VCF file.", "red")
        return

    run_snp_accession_retrieval(gene_dir, output_file)

def run_snp_accession_retrieval(gene_dir, filtered_file):
    # Set up SNPs directory and move filtered file
    live_type('Obtaining SNP accession numbers and saving to a new file.', 'magenta')
    snps_dir = os.path.join(gene_dir, 'snps')
    os.makedirs(snps_dir, exist_ok=True)
    shutil.copy(filtered_file, snps_dir)
    os.chdir(snps_dir)
    live_type(f"Copied {os.path.basename(filtered_file)} to {snps_dir}", "blue")

    # Process SNPs from filtered positions file
    new_file_name = 'processed_' + os.path.basename(filtered_file)
    process_snp_data(os.path.basename(filtered_file), new_file_name)

def process_snp_data(input_file, output_file):
    df = pd.read_csv(input_file)
    if df.empty:
        live_type(f"The DataFrame for {input_file} is empty. Skipping.", "yellow")
        return

    df.dropna(subset=['Chromosome', 'Position'], inplace=True)

    # Ensure Chromosome values are treated as strings
    df['Chromosome'] = df['Chromosome'].astype(str)
    df['Position'] = df['Position'].astype(int)

    if 'SNP' not in df.columns:
        df['SNP'] = 'N/A'

    live_type('Fetching SNP accession numbers from Ensembl...', 'magenta')
    
    for index, row in df.iterrows():
        live_type(f"Processing row {index + 1}/{len(df)}", "yellow")
        endpoint = f"http://rest.ensembl.org/overlap/region/human/{row['Chromosome']}:{row['Position']}-{row['Position']}?feature=variation"
        response = requests.get(endpoint, headers={"Content-Type": "application/json"})

        if response.status_code == 200 and response.json():
            snp_id = response.json()[0]['id']
            df.at[index, 'SNP'] = snp_id
            live_type(f"Position: {row['Position']}, SNP: {snp_id}", "green")
        else:
            live_type(f"Position: {row['Position']}, SNP: N/A", "red")

        # Save progress every 100 rows to avoid data loss
        if index % 100 == 0:
            df.to_csv(output_file, index=False)

    df.to_csv(output_file, index=False)
    live_type(f'SNP accession numbers saved as {output_file}.', 'green')

def find_csv_file():
    pattern = re.compile(r"processed_.*\.csv")
    for filename in os.listdir('.'):
        if pattern.match(filename):
            return filename
    return None

def process_csv_file(csv_file):
    df = pd.read_csv(csv_file)
    if 'SNP' not in df.columns:
        live_type(f"Error: 'SNP' column not found in {csv_file}", "red")
        return

    df_filtered = df[(df['SNP'] != 'N/A') & (df['SNP'].notna())]
    snp_list = df_filtered['SNP'].tolist()
    output_txt_file = csv_file.replace('.csv', '_snp_list.txt')
    with open(output_txt_file, 'w') as f:
        for snp in snp_list:
            f.write(f"dbsnp {snp}\n")

    live_type(f"Processed SNP list saved to {output_txt_file}", "green")

def final_snp_processing():
    csv_file = find_csv_file()
    if csv_file:
        live_type(f"Found CSV file: {csv_file}", "blue")
        process_csv_file(csv_file)
    else:
        # Not printing the "No CSV file found" message as previously requested
        pass

if __name__ == "__main__":
    filter_positions()
    final_snp_processing()
