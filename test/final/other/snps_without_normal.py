import os
import csv
import shutil
from Bio import Entrez
import pandas as pd

# Get the input file names and chromosome
merged_file = input('>>> Enter the path to the ch#_cancer_comb.csv file: ')
normal = input('>>> Enter the path to the ch#_norm_comb.csv file: ')
ch_number = input('>>> Enter the chromosome number you are analyzing 1-22, X, Y?')

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
print(f"The final merged file for chromosome {ch_number} has been saved to {output_file}")

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
        writer.writerow(row + [combined_column])
        print(combined_column)
print("\033[1;45mColumns combined and output written to\033[0m 'output.csv'")
input('Press Enter to continue...')
print("\033[1;45mQuerying NCBI for SNP accession numbers...\033[0m")

# Read the combined column from the output CSV file into a list
with open('output.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    rows = list(reader)
    col1 = [row[3] for row in rows]

# Print the list of combined column values
print(col1)

# Check if the list has a chromosome position that is not in the format chr#:pos
for position in col1:
    if ':' not in position:
        print("\033[1;45mError: The following chromosome position is not in the correct format:\033[0m", position)
        print("\033[1;45mPlease check the input file and try again.\033[0m")
        exit()

# Fetch the SNP accession numbers from NCBI
Entrez.email = "your.email@example.com"
handle = Entrez.esearch(db="snp", term=" ".join(col1))
record = Entrez.read(handle)
id_list = record["IdList"]

# Write the SNP accession numbers to a new column in the output CSV file
with open('output.csv', 'w', newline='') as csvfile_out:
    writer = csv.writer(csvfile_out)
    for i, row in enumerate(rows):
        if i == 0:
            writer.writerow(row + ['SNP Accession'])
        else:
            writer.writerow(row + [id_list[i-1]])

# Print a message to indicate that the operation is complete
print("\033[1;45mSNP accession numbers fetched and written to a new column in\033[0m 'output.csv'")
