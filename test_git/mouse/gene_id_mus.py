import requests
import json

# Get user input of a list of chromosome 11 positions
positions = input("Please enter a list of chromosome 11 positions separated by commas: ").split(",")

# Trim whitespace and convert positions to integers
positions = [int(position.strip()) for position in positions]

# Ensembl REST API URL
ensembl_api_url = "https://rest.ensembl.org"

# Function to get the gene at the specified position
def get_gene(chromosome, position):
    endpoint = f"/overlap/region/mus_musculus/{chromosome}:{position}-{position}?feature=gene;content-type=application/json"
    headers = {"Content-Type": "application/json"}

    response = requests.get(ensembl_api_url + endpoint, headers=headers)

    if response.status_code == 200:
        gene_data = json.loads(response.text)
        if len(gene_data) > 0:
            return gene_data[0]["external_name"], gene_data[0]["id"]
        else:
            return "No gene found", ""
    else:
        return "Error", ""

# Loop through the positions and find genes
for position in positions:
    gene_name, gene_id = get_gene("11", position)
    print(f"Position {position}: {gene_name} ({gene_id})")

