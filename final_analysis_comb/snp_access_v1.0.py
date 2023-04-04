#!/usr/bin/env python3

import csv
from Bio import Entrez
import shutil



# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'


# Get the name of the input CSV file from the user
csv_file = input(f"{MAGENTA}\nPlease enter the name of the CSV file:{RESET} ")

# Set the names of the source and destination files
source_file = csv_file
destination_file = 'copy_' + csv_file

# Copy the source file to the destination file
shutil.copyfile(source_file, destination_file)

# Print a message to indicate that the operation is complete
print(f"{MAGENTA}\nCopy of CSV file created{RESET}")

# Set email to identify yourself to NCBI
your_email = input(f"{MAGENTA}\nTo access NCBI, please enter your email address:{RESET} ")
Entrez.email = your_email

# Open the input and output CSV files
with open(csv_file, 'r') as csvfile_in, open('snp_accession_output.csv', 'w', newline='') as csvfile_out:
    # Create CSV reader and writer objects
    reader = csv.reader(csvfile_in)
    writer = csv.writer(csvfile_out)

    # Write the header row to the output file with a new column for SNP Accession
    header_row = next(reader)
    header_row.append('SNP Accession')
    writer.writerow(header_row)

    # Loop through the rows in the input file and add SNP accessions to a new column
    for row in reader:
        combined_column = row[0] + ':' + row[1]
        print(combined_column)
        handle = Entrez.esearch(db="snp", term=combined_column)
        record = Entrez.read(handle)
        handle.close()

        # Append the SNP accession numbers to the row
        rs_accessions = ['rs' + accession for accession in record["IdList"]]
        row.append(','.join(rs_accessions))
        print(rs_accessions)

        # Write the updated row to the output file
        writer.writerow(row)

# Print a message to indicate that the operation is complete
print(f"{MAGENTA}\nSNP Accessions added to new column in 'snp_accession_output.csv'{RESET}")



