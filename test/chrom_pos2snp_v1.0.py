#!/usr/bin/env python3

import csv
from Bio import Entrez
import shutil
import os

# Get the name of the input CSV file from the user
csv_file = input("\033[1;45mPlease enter the name of the CSV file:\033[0m")

# Set the names of the source and destination files
source_file = csv_file
destination_file = 'copy_' + csv_file

# Copy the source file to the destination file
shutil.copyfile(source_file, destination_file)

# Print a message to indicate that the operation is complete
print("Copy of CSV file created")

# Combine columns in the input CSV file and write the output to a new CSV file
with open(csv_file, 'r') as csvfile_in, open('output.csv', 'w', newline='') as csvfile_out:
    reader = csv.reader(csvfile_in)
    writer = csv.writer(csvfile_out)
    for row in reader:
        combined_column = row[1] + ':' + row[2]
        writer.writerow([combined_column])
print("\033[1;45mColumns combined and output written to\033[0m 'output.csv'")
input('Press Enter to continue...')
print("\033[1;45mQuerying NCBI for SNP accession numbers...\033[0m")

# Read the combined column from the output CSV file into a list
with open('output.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    col1 = [row[0] for row in reader]

# Print the list of combined column values
#print(col1)

# Set email to identify yourself to NCBI
your_email = input("Please enter your email address: ")
Entrez.email = your_email

#SNP accession numbers found
print('\n')
print("\033[1;45mSNP accession numbers found:\033[0m")

# Create a list to store the SNP accession numbers
rs_accessions = []

# Loop through chromosome positions and query NCBI for SNP accession numbers
for position in col1:
    handle = Entrez.esearch(db="snp", term=position)
    record = Entrez.read(handle)
    handle.close()



    # Append the SNP accession numbers to the list
    for accession in record["IdList"]:
        rs_accession = 'rs' + accession
        rs_accessions.append(rs_accession)
        print(rs_accession)

print("\033[1;45mWriting SNP accessions to CSV file...\033[0m")
# Open the input and output CSV files
with open(csv_file, 'r') as csvfile_in, open('snp_' + csv_file, 'w', newline='') as csvfile_out:
    # Create CSV reader and writer objects
    reader = csv.reader(csvfile_in)
    writer = csv.writer(csvfile_out)

    # Write the header row to the output file with a new column for SNP Accession
    header_row = next(reader)
    header_row.append('SNP Accession')
    writer.writerow(header_row)

    # Loop through the rows in the input file and add SNP accessions to a new column
    for row in reader:
        combined_column = row[1] + ':' + row[2]
        handle = Entrez.esearch(db="snp", term=combined_column)
        record = Entrez.read(handle)
        handle.close()

        # Append the SNP accession numbers to the row
        rs_accessions = ['rs' + accession for accession in record["IdList"]]
        row.append(','.join(rs_accessions))

        # Write the updated row to the output file
        writer.writerow(row)
os.remove('output.csv')

# Print a message to indicate that the operation is complete
print('\n')
print("\033[1;45mSNP Accessions added to new column in\033[0m 'snp_" + csv_file + "'")
