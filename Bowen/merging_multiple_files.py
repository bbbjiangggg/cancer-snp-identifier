#!/usr/bin/env python3

import os
import shutil
import pandas as pd
from decimal import *
import csv
from Bio import Entrez
from functools import reduce


# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# The program combines CSV files for different cancers and finds the common SNPs

# Get a list of CSV files in the current directory
data_frames = []
csvs = []
print(f'{MAGENTA}**********************************************************************{RESET}')
print(f'{MAGENTA}*****  Here are the CSV files in your current directory: {RESET}')
for i, filename in enumerate(os.listdir()):
    if filename.endswith('.csv'):
        csvs.append(filename)
        print(f'{MAGENTA}{i+1}. {filename}{RESET}')
if not csvs:
    print(f'{MAGENTA}There are no CSV files in the current directory.{RESET}')
    exit()
print(f'{MAGENTA}***** You have {len(csvs)} CSV files in your current directory.{RESET}')

# Ask user how many CSV files to merge
while True:
    repeats = input(f'{MAGENTA}>>> How many CSV files do you wish to combine (Enter a number or "ALL" to merge all): {RESET}')
    if repeats.lower() == 'all':
        repeats = len(csvs)
        break
    try:
        repeats = int(repeats)
        if repeats > 0 and repeats <= len(csvs):
            break
        else:
            print(f'{MAGENTA}Invalid input. Please enter a number between 1 and {len(csvs)} or "ALL".{RESET}')
    except ValueError:
        print(f'{MAGENTA}Invalid input. Please enter a number or "ALL".{RESET}')

# Ask user which chromosome to analyze
chromosome = input(f'{MAGENTA}>>> What is the chromosome that you are analyzing (Enter like ch1 or ch10): {RESET}')

# Merge the selected CSV files
for repeat in range(repeats):
    while True:
        try:
            file_index = int(input(f'{MAGENTA}> Enter the {repeat+1} CSV file number: {RESET}')) - 1
            if file_index >= 0 and file_index < len(csvs):
                break
            else:
                print(f'{MAGENTA}Invalid input. Please enter a number between 1 and {len(csvs)}.')
        except ValueError:
            print(f'{MAGENTA}Invalid input. Please enter a number.{RESET}')
    file = csvs[file_index]
    filename = file.replace('.csv', '').replace(chromosome, '').replace('comb', '').replace('_', '')
    df = pd.read_csv(file, index_col=0)
    df.pop('Individual')
    df.rename(columns={
        'Chrom': f'Chrom_{filename}',
        'Ref': f'Ref_{filename}',
        'Mut': f'Mut_{filename}',
        'Count': f'Count_{filename}',
    }, inplace=True)
    data_frames.append(df)

df_merged = reduce(lambda left, right: pd.merge(left, right, on=['Pos'], how='outer'), data_frames)
df_merged.replace('', float('NaN'), inplace=True)
df_merged.dropna(inplace=True)

# Ask user for the output
filename = input(f'{MAGENTA}>>> What do you want to name your merged file (add .csv at the end)? {RESET}')

#Write the merged data to the output file
with open(filename, 'w') as f:
    df_merged.to_csv(f)
print(f'{MAGENTA}Merging completed. File saved as {filename}.{RESET}')

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
