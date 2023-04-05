import csv
import requests

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Define the Ensembl URL
ensembl_url = "https://rest.ensembl.org/vep/human/id/"

# Ask the user to enter the input CSV file name
input_file_name = input(f"{MAGENTA}Enter the SNP accession CSV file's name:{RESET} ")

# Create the output CSV file name
output_file_name = "snp_effect_output.csv"

# Open the input and output CSV files
with open(input_file_name, "r") as input_file, open(output_file_name, "w", newline='') as output_file:
    reader = csv.DictReader(input_file)
    # Add a new column for the SNP effect to the output CSV file
    fieldnames = reader.fieldnames + ['SNP Effect']
    writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    writer.writeheader()
    # Process each row in the input CSV file
    for row in reader:
        # Get the rs accession number from the 'SNP Accession' column
        snp_id = row['SNP Accession'] if row['SNP Accession'] else 'N/A'
        # Construct the API URL
        url = ensembl_url + snp_id + "?content-type=application/json"
        # Open the URL and retrieve the data
        response = requests.get(url, headers={"Content-Type": "application/json"})
        # Check if the response is empty or not in JSON format
        if not response.content or not response.headers['content-type'].startswith('application/json'):
            effect_str = 'Unknown'
        else:
            # Decode the data to a JSON object
            try:
                data = response.json()
            except ValueError:
                effect_str = 'Unknown'
            # Check if the data list is empty
            if not data:
                effect_str = 'Unknown'
            else:
                # Extract the effect string using the JSON object
                if isinstance(data, list) and data[0]:
                    effect_str = data[0].get("most_severe_consequence", 'Unknown')
                else:
                    effect_str = 'Unknown'
        # Print the SNP effect
        print(f"The genetic effect of SNP {snp_id} is: {effect_str}")
        # Add the SNP effect to the row and write it to the output CSV file
        row['SNP Effect'] = effect_str
        writer.writerow(row)
print(f'{MAGENTA}SNP effect data written to{RESET} "{output_file_name}"')
