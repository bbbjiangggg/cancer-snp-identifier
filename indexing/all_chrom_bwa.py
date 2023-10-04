#!/usr/bin/env python3

import os

# BWA is used to index the reference chromosome

def download_and_index_chromosome(number):
    # Make necessary directory
    os.system('mkdir ' + number + '_bwa_ind')
    os.chdir(number + '_bwa_ind')
        
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


def start():
    chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
    for chrom in chromosomes:
        download_and_index_chromosome(chrom)


if __name__ == "__main__":
    start()
