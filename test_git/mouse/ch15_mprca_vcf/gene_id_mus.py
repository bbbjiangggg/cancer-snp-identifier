import requests
import json
import csv

# Get user input for the CSV file name
csv_file_name = input("Please enter the CSV file name: ")

# Read the data from the input CSV file
data = []
with open(csv_file_name, newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        data.append(row)

# Get chromosome number from user
valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
while True:
    chr = input('Enter the chromosome number of this analysis (1-22, X, or Y): ')
    if chr in valid_chromosomes:
        break
    else:
        print('Error: Invalid chromosome number. Please enter a valid chromosome number (1-22, X, or Y).')

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

# Add gene information to the data and print live results
for row in data:
    position = int(row["Pos"])
    gene_name, gene_id = get_gene(chr, position)
    if gene_name != "none":
        gene_function = get_gene_function(gene_id)
        row["Gene Name"] = gene_name
        row["Gene ID"] = gene_id
        row["Gene Function"] = gene_function
        print(f"Position {position}: {gene_name} ({gene_id}) - {gene_function}")
    else:
        row["Gene Name"] = "No gene found"
        row["Gene ID"] = ""
        row["Gene Function"] = ""
        print(f"Position {position}: No gene found")

# Save the updated data to a new CSV file
with open("updated_gene_id.csv", "w", newline="") as csvfile:
    fieldnames = list(data[0].keys())
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for row in data:
        writer.writerow(row)

print("Results saved to 'updated_gene_id.csv'")
