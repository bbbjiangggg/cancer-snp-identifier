#!/usr/bin/env python3

import os
from os import listdir
import shutil
import re
import requests
import io

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Set the parent directory path
parent_directory_path = os.getcwd()
print(f"{MAGENTA}\nParent directory path:{RESET} {parent_directory_path} ")

print(f"{MAGENTA}\nChecking subdirectories for{RESET} var.-final.vcf...")

# Get a list of all the subdirectories in the parent directory path
subdirectories = [d for d in os.listdir(parent_directory_path) if os.path.isdir(os.path.join(parent_directory_path, d))]

# Loop through each subdirectory
for subdir in subdirectories:
    # Check if the subdirectory starts with SRR or ERR
    if subdir.startswith("SRR") or subdir.startswith("ERR"):
        # Get the full path to the subdirectory
        subdirectory_path = os.path.join(parent_directory_path, subdir)
        # Get a list of all the files in the subdirectory
        files = os.listdir(subdirectory_path)
        # Check if any of the files end with "var.-final.vcf"
        if any(file.endswith("var.-final.vcf") for file in files):
            # If so, continue
            print(f"{MAGENTA}\nVerified:{RESET} {subdir}")
            continue
        else:
            # If not, print a message saying so
            subdirectory_name = os.path.basename(subdirectory_path)
            print(f"{MAGENTA}\nNo file ending with 'var.-final.vcf' found in subdirectory:{RESET} {subdirectory_name}")
            print(f"{MAGENTA}\nIt needs to be removed or reanalyzed.{RESET}")
            # Give option to remove
            remove = input(f"{BLUE}\nRemove? (yes/no){RESET} ")
            if remove.lower() == "yes" or remove.lower() == "y":
                # Remove the subdirectory and its contents
                shutil.rmtree(subdirectory_path)
                print(f"{MAGENTA}\nRemoved subdirectory{RESET} {subdirectory_path}.")
            else:
                print(f"{MAGENTA}\nSubdirectory not removed.{RESET}")

##########################################################################################

# This program will copy all vcf files from all SRA directories to a new directory
# Valid chromosome numbers
valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

# Get chromosome number from user
while True:
    chr = input(f'{MAGENTA}\nEnter the chromosome number of this analysis (1-22, X, or Y):{RESET} ')
    if chr in valid_chromosomes:
        break
    else:
        print(f'{RED}Error: Invalid chromosome number. Please enter a valid chromosome number (1-22, X, or Y).{RESET}')

cancer = input(f'{MAGENTA}\nEnter an abbreviation for the cancer type (e.g. pnca):{RESET} ')

directory = 'ch' + chr + '_' + cancer + '_vcf'
homedir = parent_directory_path
isecdir = os.path.join(homedir,directory)

# Make new directory
if directory in os.listdir():
    print(f'{RED}\nDirectory already exists, removing directory.{RESET}')
    shutil.rmtree(directory)
    os.mkdir(directory)
else:
    os.mkdir(directory)
    print(f'{MAGENTA}\nNew directory has been created:{RESET} ' + directory)

input(f'{BLUE}\nPress enter to continue...{RESET}')

os.chdir(isecdir)

# Copy all VCF files corresponding to the specified chromosome to the new directory
print('Copying all VCF files corresponding to chromosome ' + chr + ' to: ' + directory + '\n')

for item in listdir(homedir):
    if 'RR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            match = re.match(r"mapped_([0-9XY]+)\.var\.-final\.vcf", content)
            if match:
                file_chr = match.group(1)
                print(f"Found file: {content}, Chromosome in file: {file_chr}")
                if file_chr == chr:
                    print(f"Moving file {content} to {isecdir}")
                    shutil.move(os.path.join(subdir, content), os.path.join(isecdir, content))
                    print(f'Copied {content} corresponding to chromosome {chr} to: {directory}')
            else:
                print(f"No match for file: {content}")
    else:
        print(f"Skipping item: {item}")


print('All VCF files corresponding to chromosome ' + chr + ' have been copied to: ' + directory)
print('\n')
print(f'{MAGENTA}As backup, a copy of{RESET} ' + directory + f'{MAGENTA} is being created...{RESET}')
print('\n')

# Copy the directory as a backup
copyName = 'copy_' + directory
if copyName in os.listdir():
    print(f'{RED}\nDirectory already exists, removing directory.{RESET}')
    shutil.rmtree(copyName)
    shutil.copytree(isecdir, copyName)
else:
    shutil.copytree(isecdir, copyName)

print('A copy of ' + directory + ' has been created.')
print('\n')
input(f'{BLUE}\nPress enter to continue...{RESET}')

# Combining vcf reports
print(f'{MAGENTA} Combining all vcf files in:{RESET} ' + directory + '\n')
# Compress all .vcf files
print('Compressing all vcf files to gzip...')
os.system('ls *.vcf | xargs -n1 -P0 bgzip')
print('\n')
print('Organizing all compressed vcfgz files...')
print('\n')

# Create a directory for compressed files and move them there
os.system('mkdir isec_vcfgz_files')
os.system('mv *.gz* isec_vcfgz_files')

# Change directory to isec_vcfgz_files
os.chdir('isec_vcfgz_files')

# Create a text file with all file names for isec
os.system('ls -1 | paste -sd " " ->> isec_tools_commands.txt')

# Remove unnecessary strings from the text file
infile = "isec_tools_commands.txt"
outfile = "isec_tools_commands2.txt"
delete_string = ["isec_tools_commands.txt"]
fin = open(infile)
fout = open(outfile, "w+")
for line in fin:
    for word in delete_string:
        line = line.replace(word, "")
    fout.write(line)
fin.close()
fout.close()

# Read the existing text from file in READ mode 
with open('isec_tools_commands2.txt','r') as src:
    fline='bcftools isec -n +2 '
    # Prepending string
    oline=src.readlines()
    # Prepend the string we want to, on first line 
    oline.insert(0,fline) 
# Open the file in WRITE mode  
with open('isec_tools_commands2.txt','w') as src: 
   src.writelines(oline) 

# Create a name for the combined vcf file
vcf_file = 'ch' + chr + '_' + cancer + '_comb'
print(f'{MAGENTA}This is your combined vcf file name: {RESET}' + vcf_file)
print('\n')
print('Combining all vcf files...')
print('\n')

# Append the command to create the combined vcf file
with open('isec_tools_commands2.txt', 'a+') as file_object:
    # Append "| bgzip -c >bgzip -c > vcf_file.vcf.gz" at the end of file
    file_object.write('| bgzip -c >bgzip -c > ' + vcf_file + '.vcf.gz')

# Remove all white spaces > 1 and save it as "isec_tools_commands2.txt"
with open('isec_tools_commands2.txt', 'r') as file_object, open ('isec_tools_commands3.txt', 'w') as file_object2:
    for line in file_object:
        file_object2.write(re.sub('\s+',' ',line))

print('\n')
print('Indexing vcf files using tabix...')
# Tabix all vcfgz files
os.system('for f in ./*.vcf.gz; do tabix -p vcf -f $f;done')
print('\n')
print('Creating intersections, unions and complements (isec)...')
# Run isec on isec_tools_commands2.txt files
os.system('cat isec_tools_commands3.txt | bash')

# Unzip gzip file
os.system('gunzip ' + vcf_file + '.vcf.gz')

# Open unzipped file
print('\n')
print(f'{MAGENTA}Done! Your final combined vcf file is located here:{RESET}' + directory + '/' + vcf_file + '.vcf ') 
print('\n')
print('Open the file using the command: cat ' + directory + '/' + vcf_file + '.vcf | less')
print('\n')

os.system('rm bgzip isec_tools_commands.txt isec_tools_commands2.txt isec_tools_commands3.txt')

# Move comb.vcf file to parent directory
parent_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

for filename in os.listdir("."):
    if filename.endswith("comb.vcf"):
        # Move the file to the parent directory
        shutil.move(filename, parent_path)

# Create a dictionary of Google Drive links
link_list = {
    "1": "https://drive.google.com/file/d/1-09an3LHEJYnf-0s4S70ATis81dZPV1G/view?usp=sharing",
    "2": "https://drive.google.com/file/d/1xP8ELyW5Kr7fay0NulbgiPQOAx5Eq9E3/view?usp=sharing",
    "3": "https://drive.google.com/file/d/1-JNXwA4VvwcSOoulcmRmkLq8qYueltyQ/view?usp=sharing",
    "4": "https://drive.google.com/file/d/1n2T5xzMOqOobm5w97uthlvW5d7e3nV9M/view?usp=sharing",
    "5": "https://drive.google.com/file/d/1s8OPif3dKAesCgoOTIQGdAr9CdcUhe2q/view?usp=sharing",
    "6": "https://drive.google.com/file/d/1IzA-rUe8ELPyFtCePiCQgQeWDnOEPJxN/view?usp=sharing",
    "7": "https://drive.google.com/file/d/16eGfwwe1yjvyY6TO-hWj1Rx7YMvXuSbO/view?usp=sharing",
    "8": "https://drive.google.com/file/d/1uSqwMSvQpyBxPwiHhgUcNVkL7YdkcOHg/view?usp=sharing",
    "9": "https://drive.google.com/file/d/1chh47ez6Vqe1W0bdoKqxPHO8iOgMMiJe/view?usp=sharing",
    "10": "https://drive.google.com/file/d/12NXqBfLFJZwifp7_LLAzZFsfpE60czgN/view?usp=sharing",
    "11": "https://drive.google.com/file/d/1nsicoTeVg3t4AL592QGLa0QCiOmrhZz4/view?usp=sharing",
    "12": "https://drive.google.com/file/d/1hhX-NaqjzkpsA9cvFMj4F-1cEPCl6JmK/view?usp=sharing",
    "13": "https://drive.google.com/file/d/112aOwZiP4kOJhrBsBockjM0q2i1a1lx-/view?usp=sharing",
    "14": "https://drive.google.com/file/d/1vS_Fn_uVQvCATA402LC22QMnd5gi95pP/view?usp=sharing",
    "15": "https://drive.google.com/file/d/1-umaSF-rAFqqSwVeQI6djyW7QKxBm5A7/view?usp=sharing",
    "16": "https://drive.google.com/file/d/1cOFSi6ujbJV_I144je6PnBUmiAZL1-J0/view?usp=share_link",
    "17": "https://drive.google.com/file/d/1yhKZrhhw2YlBMySWa4L-DYaJ4UVYH-y8/view?usp=share_link",
    "18": "https://drive.google.com/file/d/1yvj7_SQE93sJ2kwceFr5ikYsFstOZ36I/view?usp=sharing",
    "19": "https://drive.google.com/file/d/1hRKEmwtw7uBkx8n_5wBWvdY2WqQjuOZz/view?usp=sharing",
    "20": "https://drive.google.com/file/d/15vb8C0tkFCczQ4uq-VWnsIQBIcxgVl9h/view?usp=sharing",
    "21": "https://drive.google.com/file/d/1kV21OWPOQdG_0mlgAIua1bJNhpf9UleN/view?usp=sharing",
    "22": "https://drive.google.com/file/d/1zqOFpGV8HEVBo_uDIqVi0_YVQOCL6LZZ/view?usp=share_link",
    "X": "https://drive.google.com/file/d/1y9ZlrPOMFz6bUctxe93x7YBchYBkqvIL/view?usp=sharing",
    "Y": "https://drive.google.com/file/d/1sNGm-4emFDuNOKnGjWozTtvLoIMXUpUg/view?usp=sharing",
    
}

while True:
    normal_file = input(f'{MAGENTA}Would you like to download a normal combined VCF file (yes or no)?{RESET} ')
    print('\n')
    if normal_file.lower() == 'yes':
        # Get chromosome number from user input
        chromosome = input(f'{MAGENTA}\nEnter the chromosome number of interest (1-22, X or Y):{RESET} ')

        # Check if chromosome number is in link list
        if chromosome in link_list:
            # Extract the file ID from the Google Drive link
            file_id = link_list[chromosome].split("/")[-2]

            # Download the VCF file as bytes
            print("Downloading chromosome " + chromosome + " normal combined VCF file...")
            url = f"https://drive.google.com/uc?id={file_id}"
            response = requests.get(url)
            vcf_bytes = io.BytesIO(response.content).read()

            # Write the VCF bytes to a file
            filename = f"ch{chromosome}_norm_comb.vcf"
            with open(filename, "wb") as f:
                f.write(vcf_bytes)
            print(f"{filename} {MAGENTA}downloaded successfully!{RESET}")
            shutil.move(filename, parent_path)
            break
        else:
            print(f"{RED}Chromosome number not found in link list.{RESET}")
    elif normal_file.lower() == 'no':
        break
    else:
        print(f"{RED}Invalid input. Please enter 'yes' or 'no'.{RESET}")
