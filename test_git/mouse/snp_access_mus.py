import csv
from Bio import Entrez

# Set email to identify yourself to NCBI
your_email = input("To access NCBI, please enter your email address: ")
Entrez.email = your_email

# Get user input of a list of chromosome 11 positions
positions = input("Please enter a list of chromosome 11 positions separated by commas: ").split(",")

# Trim whitespace and convert positions to integers
positions = [int(position.strip()) for position in positions]

# Mus musculus FVB reference sequence assembly accession number
assembly_accession = "GCF_001632505.1"

# Prepare output CSV file
with open("snp_accession_output.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Position", "SNP Accession"])

    # Loop through the positions and find SNP accessions
    for position in positions:
        query = f"{assembly_accession}[Assembly] AND 11[Chromosome] AND {position}[Base Position] AND Mus musculus[Organism]"
        handle = Entrez.esearch(db="snp", term=query)
        record = Entrez.read(handle)
        handle.close()

        # Append the SNP accession numbers
        rs_accessions = ["rs" + accession for accession in record["IdList"]]

        # Write the position and SNP accessions to the output file
        writer.writerow([position, ", ".join(rs_accessions)])

print("SNP Accessions added to 'snp_accession_output.csv'")

