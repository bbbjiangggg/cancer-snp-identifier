#!/usr/bin/env python3

import os

# BWA is used to index the reference chromosome

def create_bwa_index(chromosome_number):
    # Make necessary directory
    os.system(f"mkdir {chromosome_number}_bwa_ind")
    os.chdir(f"{chromosome_number}_bwa_ind")
        
    # Download chromosome sequence
    print(f"\033[1;45mDownloading Chromosome {chromosome_number}\033[0;0;0m")
    os.system(f"wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{chromosome_number}.fa.gz")
              
    # Indexing chromosome file
    print("\033[1;45mIndexing with BWA...\033[0;0;0m")
    os.system("bwa index *.fa.gz")
    os.system("gunzip *.fa.gz")
    print("\033[1;45mIndexing with BWA complete\033[0;0;0m")
    print("\n")
    print(f"\033[1;45mYour path to reference chromosome {chromosome_number} is:\033[0;0;0m")
    os.system("realpath *.fa")
    
def start():
    # Ask for chromosome number
    chromosome_number = input("\033[1;45mPlease enter a valid chromosome number:\033[0;0;0m ")
    print("\n")

    # Validate input
    if not chromosome_number.isdigit():
        print("\033[0;101mError: Invalid input. Please enter a number.\033[0;0;0m")
        start()
    elif int(chromosome_number) <= 0 or int(chromosome_number) > 22:
        print("\033[0;101mError: Invalid input. Please enter a number between 1 and 22.\033[0;0;0m")
        start()
    elif chromosome_number.startswith("0"):
        print("\033[0;101mError: Please do not include a '0' in front of the number.\033[0;0;0m")
        start()
    else:
        create_bwa_index(chromosome_number)

start()
