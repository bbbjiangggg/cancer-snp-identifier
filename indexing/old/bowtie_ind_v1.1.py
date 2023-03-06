#!/usr/bin/env python3

import os

# bowtie is used to index the reference chromosome

def create_directory(number):
    directory_name = f"{number}_bowtie_ind"
    if not os.path.exists(directory_name):
        os.mkdir(directory_name)
    os.chdir(directory_name)

def download_chromosome_file(number):
    print(f"\033[1;45m Downloading Chromosome {number}\033[0;0;0m")
    os.system(f"wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{number}.fa.gz")

def index_chromosome(number):
    print("\033[1;45m Indexing with Bowtie...\033[0;0;0m")
    os.system("bowtie2-build *.gz* bowtie")
    print("\033[1;45m Indexing with Bowtie complete\033[0;0;0m")
    print("\n")
    print(f"\033[1;45m Your path to the Bowtie index of chromosome {number} is:\033[0;0;0m")
    os.system("realpath bowtie")
    exit()

def start():
    number = input("\033[1;45m Please enter a valid chromosome number:\033[0;0;0m ")
    
    try:
        number = int(number)
        if 1 <= number <= 22 or number == "X" or number == "Y":
            create_directory(number)
            download_chromosome_file(number)
            index_chromosome(number)
        else:
            raise ValueError
    except ValueError:
        print("\033[0;101m NOTE: Please enter a valid chromosome number.\033[0;0;0m")
        print("\033[0;101m The number must be between 1 and 22, or 'X' or 'Y'.\033[0;0;0m")
        print("\n")
        start()

start()
