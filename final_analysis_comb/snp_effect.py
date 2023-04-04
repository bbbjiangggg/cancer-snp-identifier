import csv
import requests
import time

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Define the dbSNP URL
dbsnp_url = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"

# Define a function to retrieve the SNP effect
def get_snp_effect(rs_id):
    # Construct the API URL
    url = dbsnp_url + str(rs_id)
    # Open the URL and retrieve the data
    response = requests.get(url)
    # Check if the response is empty or not in JSON format
    if not response.content or not response.headers['content-type'].startswith('application/json'):
        return 'Unknown'
    # Decode the data to a JSON object
    data = response.json()
    # Extract the effect string using the JSON object
    effect_str = data["primary_snapshot_data"]["placements_with_allele"][0]["alleles"][0]["variation_allele"]["consequence"]
    # Return the effect string
    return effect_str

# Ask for input and output CSV file names
input_file_name = input(f'{MAGENTA}Enter "ch#_snp_accession_output.csv" file name:{RESET} ')
output_file_name = input(f'{MAGENTA}\nEnter a name for the SNP effect file (must end with .csv):{RESET} ')

# Open the input CSV file
with open(input_file_name, "r") as input_file:
    reader = csv.reader(input_file)
  

    # Define a function to print a spinning symbol for 5 seconds
    def print_spinner(text):
        symbols = ["|", "/", "-", "\\"]
        start_time = time.time()
        while time.time() - start_time < 5:
            for i in range(10):
                time.sleep(0.1)
                print(f'{MAGENTA}\r{text} {RESET}' + symbols[i % 4], end="")
        print(f'{MAGENTA}\n{text}...{RESET}')

    # Call the function to print the spinning symbol for "Contacting dbSNP"
    print_spinner("Contacting dbSNP")

    # Call the function to print the spinning symbol for "Predicting SNP effects"
    print_spinner("Predicting SNP effects")

    # Open the output CSV file
    with open(output_file_name, "w", newline="") as output_file:
        writer = csv.writer(output_file)
        # Process the header row
        header = next(reader)
        # Add the 'effect' column to the header
        header.append('effect')
        # Write the modified header to the output CSV file
        writer.writerow(header)
        # Process each row in the input CSV file
        for row in reader:
            # Get the rs accession number from the last column
            rs_id = row[-1]
            # Get the SNP effect from dbSNP
            effect = get_snp_effect(rs_id)
            # Append the SNP effect to the row
            row.append(effect)
            # Write the row to the output CSV file
            writer.writerow(row)
print(f'{MAGENTA}The SNP effect file has been saved as{RESET} {output_file_name}')