#!/usr/bin/env python3

import os
import shutil
import pandas as pd
from decimal import *
import csv
from Bio import Entrez


# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

print('\n')
print(f'{MAGENTA}This program will merge the cancer and normal csv files.{RESET}')
input(f'{BLUE}\nPress enter to continue...{RESET}')

#cwd = os.getcwd()
#print(cwd)

#change directory to csv_files
os.chdir("csv_files")

#These are the files in the current directory
current_directory = os.getcwd()
files = os.listdir(current_directory)

print("These are the files in the current directory:")
for file in files:
    print(file)
    

# Get the input file names and chromosome
merged_file = input(f'{MAGENTA}\nEnter the name of the ch#_cancer_comb.csv file:{RESET} ')
normal = input(f'{MAGENTA}\nEnter the name of the ch#_norm_comb.csv file:{RESET} ')
ch_number = input(f'{MAGENTA}\nEnter the chromosome number you are analyzing 1-22, X, Y?{RESET} ')

# Load the normal and merged files into pandas dataframes
df_normal = pd.read_csv(normal, index_col=[0])
df_merged = pd.read_csv(merged_file, index_col=[0])

# Get a list of positions from the normal file
positions = df_normal['Pos'].tolist()

# Filter the merged dataframe to remove normal SNPs
df_filtered = df_merged.loc[df_merged['Pos'].isin(positions)==False]

# Save the filtered dataframe to a new CSV file
output_file = 'ch' + ch_number + '_final_merge.csv'
df_filtered.to_csv(output_file, index=False)

# Print a message to indicate that the operation is complete
print(f"{MAGENTA}The final merged file for chromosome {ch_number} has been saved as {output_file}{RESET}")
print('\n')
##################################################################################

print(f'{MAGENTA}This program will access NCBI and find SNP accession numbers.{RESET}')
input(f'{BLUE}\nPress enter to continue...{RESET}')

#These are the files in the current directory
current_directory = os.getcwd()
files = os.listdir(current_directory)

print("These are the files in the current directory:")
for file in files:
    print(file)


# Get the name of the input CSV file from the user
csv_file = input(f'{MAGENTA}\nPlease enter the name of the "ch#_final_merge.csv" file:{RESET} ')

# Set the names of the source and destination files
source_file = csv_file
destination_file = 'copy_' + csv_file

# Copy the source file to the destination file
shutil.copyfile(source_file, destination_file)

# Print a message to indicate that the operation is complete
print(f"{MAGENTA}\nCopy of csv file created{RESET}")

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