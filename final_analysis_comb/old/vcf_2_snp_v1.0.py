#!/usr/bin/env python3

import os
from os import listdir
import shutil
import re
import requests
import io
import pandas as pd
from decimal import *
import numpy as np
import csv
from Bio import Entrez


# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

print(f"{MAGENTA}\nThis program will copy all vcf files from all SRA directories to a new directory{RESET}")
input(f'{BLUE}\nPress enter to continue...{RESET}')

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
                os.system(f"rm -r {subdirectory_path}")
                print(f"{MAGENTA}\nRemoved subdirectory{RESET} {subdirectory_path}.")
            else:
                print(f"{MAGENTA}\nSubdirectory not removed.{RESET}")

##########################################################################################

#This program will copy all vcf files from all SRA directories to a new directory
#chr = input(f'{MAGENTA}\nEnter the chrosome number of this analysis:{RESET} ')
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
#print(f"{MAGENTA}\nThis is your directory's name:{RESET} " + directory)

homedir = parent_directory_path


isecdir = os.path.join(homedir,directory)

#make new directory
if directory in os.listdir():
    print(f'{RED}\nDirectory already exists, removing directory.{RESET}')
    os.system('rm -r ' + directory)
    os.mkdir(directory)
else:
    os.mkdir(directory)
    print(f'{MAGENTA}\nNew directory has been created:{RESET} ' + directory)
    
input(f'{BLUE}\nPress enter to continue...{RESET}')
print('\n') 

os.chdir(isecdir)
#copy all vcf files to new directory
print('Copying all vcf files to: ' + directory + '\n')

for item in listdir(homedir):
    if 'RR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                shutil.move(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('All vcf files have been copied to: ' + directory)
print('\n')
print(f'{MAGENTA}As backup, a copy of{RESET} ' + directory + f'{MAGENTA} is being created...{RESET}')
print('\n')

copyName = 'copy_'+directory
if copyName in os.listdir():
    print(f'{RED}\nDirectory already exists, removing directory.{RESET}')
    os.system('rm -r ' + copyName)
    os.mkdir(copyName)
else:
    #copydir = os.path.join(homedir, copyName)
    shutil.copytree(isecdir, copyName)

print('A copy of ' + directory + ' has been created.')
print('\n')
input(f'{BLUE}\nPress enter to continue...{RESET}')
print('\n')

#Combining vcf reports
print(f'{MAGENTA} Combining all vcf files in:{RESET} ' + directory + '\n')
#this command will bgzip all .vcf files, this program must be located in the same directory
#if all files are already bgzipped then it will throw back:
#"ls: cannot access '*.vcf': No such file or directory"

print('Compressing all vcf files to gzip...')


os.system('ls *.vcf | xargs -n1 -P0 bgzip')
print('\n')
print('Organizing all compressed vcfgz files...')
print('\n')

#move all .vcf files to a directory called isec_vcfgz_files

os.system('mkdir isec_vcfgz_files')
os.system('mv *.gz* isec_vcfgz_files')

#change directory to isec_vcfgz_files
os.chdir('isec_vcfgz_files')



#this command will Write all file names into a txt file, same line, 
#one space, named "isec_tools_commands.txt"
os.system('ls -1 | paste -sd " " ->> isec_tools_commands.txt')

#deleting string that results from previous command
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



#read the existing text from file in READ mode 
with open('isec_tools_commands2.txt','r') as src:
    fline='bcftools isec -n +2 '
    #Prepending string
    oline=src.readlines()
    #prepend the string we want to, on first line 
    oline.insert(0,fline) 
#open the file in WRITE mode  
with open('isec_tools_commands2.txt','w') as src: 
   src.writelines(oline) 

#vcf_file = input('\033[1;45m Enter the name you wish to give the combined vcf files (e.g. ch22_pnca_comb): \033[0m')

vcf_file = 'ch' + chr + '_' + cancer + '_comb'
print(f'{MAGENTA}This is your combined vcf file name: {RESET}' + vcf_file)
print('\n')
print('Combining all vcf files...')
print('\n')

#open a file with access mode 'a+'web: https://stackabuse.com/file-handling-in-python/
with open('isec_tools_commands2.txt', 'a+') as file_object:
    # Append "| bgzip -c >bgzip -c > vcf_file.vcf.gz" at the end of file
    file_object.write('| bgzip -c >bgzip -c > ' + vcf_file + '.vcf.gz')
    
#remove all white spaces > 1 and save it as "isec_tools_commands2.txt"
with open('isec_tools_commands2.txt', 'r') as file_object, open ('isec_tools_commands3.txt', 'w') as file_object2:
    for line in file_object:
        file_object2.write(re.sub('\s+',' ',line))

print('\n')
print('Indexing vcf files using tabix...')


#tabix all vcfgz files
os.system('for f in ./*.vcf.gz; do tabix -p vcf -f $f;done')
print('\n')
print('Creating intersections, unions and complements (isec)...')

#run isec on isec_tools_commands2.txt files
os.system('cat isec_tools_commands3.txt | bash')

#unzip gzip file
os.system('gunzip '+ vcf_file + '.vcf.gz')

#open unzipped file
print('\n')
print(f'{MAGENTA}Done! Your final combined vcf file is located here:{RESET}' + directory + '/' + vcf_file + '.vcf ') 
print('\n')
print('Open the file using the command: cat ' + directory + '/' + vcf_file + '.vcf | less')
print('\n')


os.system('rm bgzip isec_tools_commands.txt isec_tools_commands2.txt isec_tools_commands3.txt')

#move comb.vcf file to parent directory
parent_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

for filename in os.listdir("."):
    if filename.endswith("comb.vcf"):
        # Move the file to the parent directory
        shutil.move(filename, parent_path)

print(f'{MAGENTA}You will now have the opportunity to download the "normal_comb.vcf"{RESET} ')
input(f'{BLUE}\nPress enter to continue...{RESET}')
print('\n')

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

#######################################################################################################
print('\n')
print(f"{MAGENTA}This program will convert VCF files to CSV files.{RESET}")
input(f'{BLUE}\nPress enter to continue...{RESET}')

os.chdir("..")  # move one level up in the directory tree
os.chdir("..")  # move one level up in the directory tree

# get the current working directory
cwd = os.getcwd()


# list all directories in the current working directory
dir_list = os.listdir(cwd)

# loop through the directories and find the first one that ends with "_vcf"
vcf_dir = None
for d in dir_list:
    if os.path.isdir(d) and d.endswith("_vcf"):
        # set vcf_dir to the first directory that ends with "_vcf"
        vcf_dir = d
        break

if vcf_dir is None:
    print(f"{RED}No directory ending with '_vcf' found in current directory.{RESET}")
else:
    # change the current working directory to the first directory that ends with "_vcf"
    os.chdir(vcf_dir)
    print(f"{MAGENTA}\nChanged current directory to:{RESET} ", os.getcwd())

print('\n')

#os.chdir(directory + '/isec_vcfgz_files/')

def snp_above_70(filename):
    getcontext().prec = 2
    count = 0
    individual_count = 0
    allcols = []

    #check the count percentage
    with open(filename) as file:
        for line in file:
            cols = line.split('\t')
            for letter in cols[4]:
                individual_count = individual_count + 1
                if letter == "1":
                    count += 1
            if Decimal(count)/Decimal(individual_count) >= 0.70:            
                cols.append(str(count))
            allcols.append(cols)
            count = 0
            individual_count = 0


        #delete rows with NaN and sort the other counts
        df = pd.DataFrame(allcols)
        df.columns = ['Chrom', 'Pos', 'Ref', 'Mut', 'Individual', 'Count']
        df.Count = pd.to_numeric(df.Count, errors='coerce')
        df['Count'].replace('', np.nan, inplace=True)
        df.dropna(subset=['Count'], inplace=True)
        sort = df.sort_values(by = 'Count', ascending=False)
        filename = filename.replace('.vcf', '')
        sort.to_csv(filename+'.csv')
        print(sort)

print(f'{MAGENTA}Here is/are the vcf file(s) in your current directory:{RESET}')
num = 0
transfer = []
for i in os.listdir():
    if '.vcf' in i:
        print(i)
        num = num + 1
        transfer.append(i)
print(f'{MAGENTA}\nThere is/are{RESET} {num} {MAGENTA}vcf file(s) in your current directory.{RESET}')

vcfrun = input(f'{MAGENTA}Enter how many "comb.vcf" files you wish to convert to csv (type "ALL" to transfer all): {RESET}')

if vcfrun == 'ALL':
    for file in transfer:
        snp_above_70(file)
        print(f'{MAGENTA} {file} converted to csv format {RESET}')


if vcfrun != 'ALL':
    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
    vcfrun = int(vcfrun)
    vcfrun+=1
    placement = [ordinal(n) for n in range(1, vcfrun)]

    vcfs = []

    for time in range(int(vcfrun)-1):
        filename = input("Enter the " + placement[time] + " vcf file you wish to convert (type q to quit): ")
        if filename == 'q':
            break
        else: vcfs.append(filename)

    for name in vcfs:   
        snp_above_70(name)
        print(f'{MAGENTA}{name} converted to csv format {RESET}')

print('\n')

homedir = os.getcwd()

if 'csv_files' not in os.listdir():
        os.makedirs(os.path.join(homedir, 'csv_files'))
        print(f'{MAGENTA}Creating "csv" directory...{RESET}')

for file in os.listdir():
    if '.csv' in file:
        dst = os.path.join(homedir, 'csv_files')
        shutil.move(os.path.join(homedir, file),
                    os.path.join(dst, file))
print(f'{MAGENTA}\nAll vcf files have been moved into the "csv" directory. {RESET}')

#######################################################################################
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