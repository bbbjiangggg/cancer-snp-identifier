#!/usr/bin/env python3

import os

# BWA is used to index the reference chromosome

def download_and_index_chromosome(number):
    directory_name = number + '_bwa_ind'
    
    # Check if the BWA index already exists for the chromosome
    if index_exists(directory_name, number):
        print(f'\033[1;45mBWA index for chromosome {number} already exists. Skipping...\033[0;0;0m')
        return

    # Make necessary directory
    if not os.path.exists(directory_name):
        os.system('mkdir ' + directory_name)

    os.chdir(directory_name)
        
    # Download chromosome sequence
    print('\033[1;45mDownloading Chromosome '+ number + '\033[0;0;0m')
    os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.'+ number +'.fa.gz')
              
    # Indexing chromosome file
    print('\033[1;45mIndexing with BWA...\033[0;0;0m')
    os.system('bwa index *.gz*')
    os.system('bgzip -d *.gz*')
    
    # Running samtools command to faidx the chromosome
    os.system('samtools faidx Homo_sapiens.GRCh38.dna.chromosome.' + number + '.fa')
    
    print()
    print('\033[1;45mYour path to reference chromosome ' + number + ' is:\033[0;0;0m')
    os.system('realpath *.fa')
    
    # Move back to parent directory for next chromosome
    os.chdir("..")


def index_exists(directory, chrom_number):
    """
    Checks if the BWA index exists for a given chromosome directory.

    Args:
        directory: The directory where the FASTA file and BWA index are located.
        chrom_number: The chromosome number.
    
    Returns:
        bool: True if the BWA index exists, False otherwise.
    """
    # BWA creates multiple index files with extensions like .amb, .ann, etc. 
    # We can check for one of these to determine if the index exists.
    return os.path.exists(os.path.join(directory, f"Homo_sapiens.GRCh38.dna.chromosome.{chrom_number}.fa.amb"))


def start():
    chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
    for chrom in chromosomes:
        download_and_index_chromosome(chrom)


if __name__ == "__main__":
    start()
