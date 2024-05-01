
import os

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def replace_text(file_path, old_text, new_text):
    with open(file_path, 'r+') as file:
        text = file.read().replace(old_text, new_text)
        file.seek(0)
        file.write(text)
        file.truncate()

def replace_in_untrimmed_bash_srr(old_text, new_text):
    replace_text('untrimmed_bash_sra_v1.2.txt', old_text, new_text)

valid_chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
chromosome = ""
while chromosome not in valid_chromosomes:
    chromosome = input(f'{MAGENTA}Enter the chromosome number or name (e.g., 1, 2, ... 22, X, Y) to be analyzed: {RESET}')
    if chromosome not in valid_chromosomes:
        print(f"{RED}Invalid chromosome number or name. Please enter a valid chromosome.{RESET}")

# Construct the paths for BWA and Bowtie files based on the chromosome
bwa_chrom_path = f"/usr/local/bin/bwa/{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
bowtie_index_path = f"/usr/local/bin/bowtie/{chromosome}_bowtie_ind/bowtie"

# Determine the path for the trimmomatic-0.39.jar file
jar_file = '~/Trimmomatic-0.39/trimmomatic-0.39.jar'
jar_path = os.path.expanduser(jar_file)

if os.path.exists(jar_path):
    trim_path = jar_path
    replace_in_untrimmed_bash_srr('trim_path', trim_path)
    print(f'{MAGENTA}3){RESET} {jar_path} is the absolute path.')
else:
    print(f'NOTE: {jar_path} does not match your absolute path.')
    print('You have a different path for trimmomatic-0.39.jar')
    trim_path = input(f'{MAGENTA}3){RESET} Copy and paste the absolute path to your trimmomatic-0.39.jar file: ')
    replace_in_untrimmed_bash_srr('trim_path', trim_path)
