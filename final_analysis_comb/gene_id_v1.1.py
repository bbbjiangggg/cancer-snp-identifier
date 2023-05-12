import csv
import requests

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Define the Ensembl URL
ensembl_url = "https://rest.ensembl.org/"

# Ask the user to enter the input CSV file name
input_file_name = input(f"{MAGENTA}Enter the SNP accession CSV file's name:{RESET} ")

# Create the output CSV file name
output_file_name = "snp_gene_output.csv"

# Open the input and output CSV files
with open(input_file_name, "r") as input_file, open(output_file_name, "w", newline='') as output_file:
    reader = csv.DictReader(input_file)
    # Add a new column for the gene name to the output CSV file
    fieldnames = reader.fieldnames + ['Gene Name']
    writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    writer.writeheader()
    # Process each row in the input CSV file
    for row in reader:
        # Get the SNP accession from the 'Accession' column
        snp_accession = row['SNP Accession'] if row['SNP Accession'] else 'N/A'
        # Construct the API URL to get the SNP position
        url = ensembl_url + "variation/human/" + snp_accession + "?content-type=application/json"
        # Open the URL and retrieve the data
        response = requests.get(url, headers={"Content-Type": "application/json"})
        # Check if the response is empty or not in JSON format
        if not response.content or not response.headers['content-type'].startswith('application/json'):
            snp_position = 'N/A'
        else:
            # Decode the data to a JSON object
            try:
                data = response.json()
            except ValueError:
                snp_position = 'N/A'
            # Check if the data contains the SNP position
            if 'mappings' not in data or not data['mappings']:
                snp_position = 'N/A'
            else:
                # Extract the SNP position using the JSON object
                snp_position = str(data['mappings'][0]['start'])
        # Construct the API URL to get the gene name
        url = ensembl_url + "overlap/region/human/" + snp_position + "?feature=gene;content-type=application/json"
        # Open the URL and retrieve the data
        response = requests.get(url, headers={"Content-Type": "application/json"})
        # Check if the response is empty or not in JSON format
        if not response.content or not response.headers['content-type'].startswith('application/json'):
            gene_name = 'Unknown'
        else:
            # Decode the data to a JSON object
            try:
                data = response.json()
            except ValueError:
                gene_name = 'Unknown'
            # Check if the data list is empty
            if not data:
                gene_name = 'Unknown'
            else:
                # Extract the gene name using the JSON object
                if isinstance(data, list) and data[0]:
                    gene_name = data[0].get("gene_name", 'Unknown')
                else:
                    gene_name = 'Unknown'
        

        # Print the gene name
        print(f"The gene at SNP position {snp_position} is: {gene_name}")
        # Add the gene name to the row and write it to the output CSV file
        row['Gene Name'] = gene_name
        writer.writerow(row)
print(f'{MAGENTA}Gene name data written to{RESET} "{output_file_name}"')

