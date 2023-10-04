#!/usr/bin/env python3

import os
import subprocess


def make_directory(number):
    """
    Creates a new directory for the Bowtie index.

    Args:
        number: A string representing the chromosome number.
    """
    directory_name = f"{number}_bowtie_ind"
    os.makedirs(directory_name, exist_ok=True)
    return directory_name


def download_chromosome_file(directory, number):
    """
    Downloads the FASTA file for the specified chromosome number.

    Args:
        directory: The directory where the FASTA file will be saved.
        number: A string representing the chromosome number.
    """
    url = f"http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{number}.fa.gz"
    subprocess.run(["wget", "-P", directory, url])


def index_chromosome(directory, number):
    """
    Creates a Bowtie index for the specified chromosome.

    Args:
        directory: The directory where the FASTA file is located.
        number: A string representing the chromosome number.
    """
    print("\033[1;45m Indexing with Bowtie...\033[0;0;0m")
    subprocess.run(["bowtie2-build", os.path.join(directory, f"Homo_sapiens.GRCh38.dna.chromosome.{number}.fa.gz"), os.path.join(directory, "bowtie")])
    
    print("\n")
    print(f"\033[1;45m Your path to the Bowtie index of chromosome {number} is:\033[0;0;0m")
    bowtie_path = subprocess.run(["realpath", os.path.join(directory, "bowtie")], capture_output=True)
    print(bowtie_path.stdout.decode('utf-8'))
    print("\n")   


def index_exists(directory):
    """
    Checks if the Bowtie index exists for a given chromosome directory.

    Args:
        directory: The directory where the FASTA file and Bowtie index are located.
    
    Returns:
        bool: True if the Bowtie index exists, False otherwise.
    """
    # Bowtie2 creates multiple index files with extensions .bt2 or .bt2l. 
    # We can check for one of these to determine if the index exists.
    return os.path.exists(os.path.join(directory, "bowtie.1.bt2"))


def index_all_chromosomes():
    """
    Downloads and creates a Bowtie index for all human chromosomes.
    """
    chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
    
    for chrom in chromosomes:
        directory = make_directory(chrom)
        
        # Check if the index already exists
        if index_exists(directory):
            print(f"\033[1;45m Bowtie index for chromosome {chrom} already exists. Skipping...\033[0;0;0m")
            continue
        
        download_chromosome_file(directory, chrom)
        index_chromosome(directory, chrom)


def main():
    """
    The main function that runs the program.
    """
    index_all_chromosomes()


if __name__ == "__main__":
    main()
