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
    os.chdir(directory_name)


def download_chromosome_file(number):
    """
    Downloads the FASTA file for the specified chromosome number.

    Args:
        number: A string representing the chromosome number.
    """
    url = f"http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{number}.fa.gz"
    subprocess.run(["wget", url])


def index_chromosome(number):
    """
    Creates a Bowtie index for the specified chromosome.

    Args:
        number: A string representing the chromosome number.
    """
    print("\033[1;45m Indexing with Bowtie...\033[0;0;0m")
    subprocess.run(["bowtie2-build", f"Homo_sapiens.GRCh38.dna.chromosome.{number}.fa.gz", "bowtie"])
    
    print("\n")
    print(f"\033[1;45m Your path to the Bowtie index of chromosome {number} is:\033[0;0;0m")
    bowtie_path = subprocess.run(["realpath", "bowtie"])
    print(bowtie_path)
    print("\n")   

    
def start():
    """
    Prompts the user for a chromosome number and creates a Bowtie index for the corresponding chromosome.
    """
    number = input("\033[1;45m Please enter a valid chromosome number (1-22, X, Y):\033[0;0;0m ")

    if number.isdigit() and int(number) in range(1, 23) or number.upper() in ["X", "Y"]:
        number = str(number)
        make_directory(number)
        download_chromosome_file(number)
        index_chromosome(number)
    else:
        print("\033[0;101m NOTE: Please enter a valid chromosome number.\033[0;0;0m")
        print("\033[0;101m The number must be between 1 and 22, or 'X' or 'Y'.\033[0;0;0m")
        print("\n")
        start()


def main():
    """
    The main function that runs the program.
    """
    start()


if __name__ == "__main__":
    main()
