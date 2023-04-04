#!/usr/bin/env python3

import os
import shutil
import pandas as pd
from decimal import *
import csv
from Bio import Entrez
from functools import reduce
import requests
import time
from termcolor import colored
import importlib
import subprocess


# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'



packages = ['os', 'shutil', 'pandas', 'decimal', 'csv', 'Bio', 'functools', 'requests', 'termcolor']

# Check if packages are installed, install them if necessary
for package in packages:
    try:
        importlib.import_module(package)
    except ImportError:
        print(f'{package} is not installed. Installing...')
        subprocess.run(['pip', 'install', package])



#The program combines CSV files for different cancers and finds the common SNPs

dataCount = 0
data_frames = []
csvs = []

num = 0
for i in os.listdir():
    if '.csv' in i:
        csvs.append(i)
        print(i)
        num += 1
print(f'{MAGENTA}You have {num} csv files in your current directory.{RESET}')

repeats = input(f'{MAGENTA}\nHow many csv files do you wish to combine (Enter a number or "ALL" to merge all): {RESET}')
chromosome = input(f'{MAGENTA}\nEnter the chromosome number you are analyzing (1-22, X or Y): {RESET}')


if repeats == 'ALL':
    for repeat in range(num):
        globals()['df%s'%(repeat+1)] = pd.read_csv(csvs[repeat], index_col = [0])
        globals()['df%s'%(repeat+1)].pop('Individual')

        file = csvs[repeat]
        
        filename = file.replace('.csv', '')
        filename = filename.replace(chromosome, '')
        filename = filename.replace('comb', '')
        filename = filename.replace('_', '')

        globals()['df%s'%(repeat+1)].rename(columns={"Chrom": ("Chrom_" + str(filename)), "Ref": ("Ref_" + str(filename)),
                                                 "Mut": ("Mut_" + str(filename)), "Count": ("Count_" + str(filename))}, inplace = True)

        dataCount = dataCount + 1
        
if repeats != 'ALL':
    for repeat in range(int(repeats)):
        file = input('\033[1;35m> Enter the '+ str(repeat+1) + ' csv file: \033[0m')

        filename = file.replace('.csv', '')
        filename = filename.replace(chromosome, '')
        filename = filename.replace('comb', '')
        filename = filename.replace('_', '')
        
        globals()['df%s'%(repeat+1)] = pd.read_csv(file, index_col = [0])
        globals()['df%s'%(repeat+1)].pop('Individual')
        
        globals()['df%s'%(repeat+1)].rename(columns={"Chrom": ("Chrom_" + str(filename)), "Ref": ("Ref_" + str(filename)),
                                                     "Mut": ("Mut_" + str(filename)), "Count": ("Count_" + str(filename))}, inplace = True)
        
        dataCount = dataCount + 1

for num in range(dataCount):
    data_frames.append(globals()['df%s'%(num+1)])

df_merged = reduce(lambda left,right: pd.merge(left, right, on=['Pos'], how='outer'), data_frames)

nan_value = float("NaN")
df_merged.replace("", nan_value, inplace=True)

for column in df_merged.columns.tolist():
    df_merged.dropna(subset = [column], inplace=True)

wish = input(f'{MAGENTA}\nEnter a name for your merged file (must end with ".csv")? {RESET}')


# Replace first column with "chr" and "chromosome" value
df_merged.insert(0, "chr", chromosome)
df_merged = df_merged.reset_index(drop=True)
df_merged.to_csv(wish)


# Open the merged CSV file and delete the first column
df = pd.read_csv(wish)
df = df.drop(columns=['Unnamed: 0'])
df.to_csv(wish, index=False)



print(f'{MAGENTA}\nMerging Completed {RESET}')


##################################################################################

print(f'{MAGENTA}\nThis program will now access NCBI and find SNP accession numbers.{RESET}')
input(f'{BLUE}\nPress enter to continue...{RESET}')

#These are the files in the current directory
current_directory = os.getcwd()
files = os.listdir(current_directory)

print("These are the files in the current directory:")
for file in files:
    print(file)


# Get the name of the input CSV file from the user
print(f'\n{colored("Is this the name of the csv file for SNP analysis (yes/no)?", "magenta")} {wish}')
if input().lower() == 'yes':
    csv_file = wish
else:
    csv_file = input(f'{colored("Please enter the name of the csv file for SNP analysis:", "magenta")} ')

# Validate the input CSV file name
while not csv_file.endswith('.csv'):
    print(colored('Input file name must end with ".csv".', 'red'))
    csv_file = input(f'{colored("Please enter the name of the csv file for SNP analysis:", "magenta")} ')


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

print(f"{MAGENTA}\nFinding SNP accession numbers...{RESET}")

# Open the input and output CSV files
with open(csv_file, 'r') as csvfile_in, open(f'ch{chromosome}_snp_accession_output.csv', 'w', newline='') as csvfile_out:
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
snp_file = f'ch{chromosome}_snp_accession_output.csv'
# Print a message to indicate that the operation is complete
print(f"{MAGENTA}\nSNP Accessions added to new column in{RESET} 'ch{chromosome}_snp_accession_output.csv'")

#####################################################################
print(f'{MAGENTA}\nThis program will now access dbSNP and find SNP effects.{RESET}')
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

# Import necessary modules
from termcolor import colored

# Ask for input and output CSV file names
print(f'\n{colored("Is this your input file name (yes/no)?", "magenta")} {snp_file}')
if input().lower() == 'yes':
    input_file_name = snp_file
else:
    input_file_name = input(f'{colored("Enter the name of the SNP accession file (must end with .csv):", "magenta")} ')

output_file_name = input(f'\n{colored("Enter a name for the SNP effect file (must end with .csv):", "magenta")} ')

# Validate input and output file names
while not input_file_name.endswith('.csv'):
    print(colored('Input file name must end with ".csv".', 'red'))
    input_file_name = input(f'{colored("Enter the name of the SNP accession file (must end with .csv):", "magenta")} ')

while not output_file_name.endswith('.csv'):
    print(colored('Output file name must end with ".csv".', 'red'))
    output_file_name = input(f'{colored("Enter a name for the SNP effect file (must end with .csv):", "magenta")} ')

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
print(f'{MAGENTA}\nThe SNP effect file has been saved as{RESET} {output_file_name}')
