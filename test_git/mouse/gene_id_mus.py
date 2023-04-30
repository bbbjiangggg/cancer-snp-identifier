import requests
import json
import csv

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
            return "none", ""
    else:
        return "Error", ""

# Function to get gene function
def get_gene_function(gene_id):
    endpoint = f"/lookup/id/{gene_id}?content-type=application/json"
    headers = {"Content-Type": "application/json"}

    response = requests.get(ensembl_api_url + endpoint, headers=headers)

    if response.status_code == 200:
        gene_info = json.loads(response.text)
        if "description" in gene_info:
            return gene_info["description"]
        else:
            return "No description available"
    else:
        return "Error"

# Open a new CSV file for writing
with open("gene_id.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Position", "Gene Name", "Gene ID", "Gene Function"])

    # Loop through the positions and find genes
    for position in positions:
        gene_name, gene_id = get_gene("11", position)
        if gene_name != "none":
            gene_function = get_gene_function(gene_id)
            writer.writerow([position, gene_name, gene_id, gene_function])
            print(f"Position {position}: {gene_name} ({gene_id}) - {gene_function}")
        else:
            writer.writerow([position, "No gene found", "", ""])
            print(f"Position {position}: No gene found")

print("Results saved to 'gene_id.csv'")
