
#!/usr/bin/env python3

import os
import shutil
import subprocess
from pathlib import Path
import signal
import importlib

# ... (previous code remains unchanged)

# Ask the user for chromosomes to analyze
chromosomes_input = input(f"{MAGENTA}Enter the chromosomes to be analyzed, separated by a comma (e.g., 1,2,X,Y): {RESET}")
chromosomes = [chrom.strip() for chrom in chromosomes_input.split(",")]

# Loop through each chromosome for analysis
for chromosome in chromosomes:
    print(f"\n{MAGENTA}Analyzing chromosome: {chromosome}{RESET}")
    
    # Construct the paths for BWA and Bowtie files based on the chromosome
    bwa_chrom_path = f"/usr/local/bin/bwa/{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
    bowtie_index_path = f"/usr/local/bin/bowtie/{chromosome}_bowtie_ind/bowtie"

    # Print the paths
    print(f"{MAGENTA}Bowtie Index Path: {RESET}{bowtie_index_path}")
    print(f"{MAGENTA}BWA Chromosome Path: {RESET}{bwa_chrom_path}")

    # Add the path to where bowtie files are found (must end in 'bowtie')
    replace_in_untrimmed_bash_srr('bowtie_index_path', bowtie_index_path)

    # Add the path to where reference chromosome is found
    replace_in_untrimmed_bash_srr('bwa_chrom_path', bwa_chrom_path)

    # ... (rest of the code remains unchanged, excluding the email sending part)

# After all chromosomes have been analyzed, send the email
print('Sending email to ' + user + ' ....')
os.system('sendemail -f sudoroot1775@outlook.com -t ' + user + ' -u ' + job_title + '_Analysis Done -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')
