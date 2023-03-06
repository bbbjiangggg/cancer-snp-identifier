#!/usr/bin/env python3

import os
import shutil

# Get input from user
chromosome_number = input('Enter the chromosome number for this analysis: ')
cancer_type = input('Enter an abbreviation for the cancer type (e.g. pnca): ')

# Define directories and paths
directory_name = f'ch{chromosome_number}_{cancer_type}_vcf'
print(f"This is your directory's name: {directory_name}")
print("Note: Use the command 'realpath name_of_file/dir' to get the complete path.")
homedir_path = input('Enter the path where the SRA directories are currently stored: ')
directory_path = os.path.join(homedir_path, directory_name)

# Copy all vcf files to new directory
if os.path.exists(directory_name):
    print(f"Directory already exists, removing directory: {directory_name}")
    shutil.rmtree(directory_name)
os.mkdir(directory_name)
print(f"New directory has been created: {directory_name}")

print("Copying all vcf files to the new directory...")
for root, _, files in os.walk(homedir_path):
    for filename in files:
        if filename.endswith('.vcf') and 'RR' in root:
            file_path = os.path.join(root, filename)
            shutil.copy(file_path, directory_path)

print(f"All vcf files have been copied to: {directory_name}")

# Create a backup copy of the directory
copy_name = f'copy_{directory_name}'
if os.path.exists(copy_name):
    print(f"Directory already exists, removing directory: {copy_name}")
    shutil.rmtree(copy_name)
shutil.copytree(directory_path, copy_name)
print(f"A backup copy of {directory_name} has been created.")

# Combine vcf files
print("Combining all vcf files...")
os.chdir(directory_path)

# Compress all vcf files to gzip
print("Compressing all vcf files to gzip...")
os.system('ls *.vcf | xargs -n1 -P0 bgzip')

# Move all .vcf.gz files to a directory called "isec_vcfgz_files"
print("Organizing all compressed vcfgz files...")
os.mkdir('isec_vcfgz_files')
os.system('mv *.gz* isec_vcfgz_files')
os.chdir('isec_vcfgz_files')

# Write all file names into a txt file, same line, one space, named "isec_tools_commands.txt"
print("Writing vcf file names into isec_tools_commands.txt...")
with open('isec_tools_commands.txt', 'w') as file:
    for root, _, files in os.walk('.'):
        for filename in files:
            if filename.endswith('.vcf.gz'):
                file.write(filename + ' ')
    file.write('\n')

# Prepend bcftools command to the txt file
print("Prepending bcftools command to isec_tools_commands.txt...")
with open('isec_tools_commands.txt', 'r+') as file:
    contents = file.read()
    file.seek(0, 0)
    file.write(f'bcftools isec -n +2 {contents}')

print("Done!")
