#!/usr/bin/env python3

import os

# BWA is used to index the reference chromosome

def start():
    # make necessary directory
    while True:
        number = input('\033[1;45mPlease enter a valid chromosome number (1-22, X, Y): \033[0;0;0m').strip().upper()
        
        if number in ['X', 'Y'] or number.isdigit() and 1 <= int(number) <= 22:
            print('\n')
            break
        else:
            print('\033[0;101mInvalid chromosome number. Please enter a valid chromosome number (1-22, X, Y).\033[0;0;0m')
            print()
            continue
        
    os.system('mkdir ' + number + '_bwa_ind')
    os.chdir(number + '_bwa_ind')
        
    # download chromosome sequence
    print('\033[1;45mDownloading Chromosome '+ number + '\033[0;0;0m')
    os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.'+ number +'.fa.gz')
              
    # indexing chromosome file
    print('\033[1;45mIndexing with BWA...\033[0;0;0m')
        
    os.system('bwa index *.gz*')
    os.system('bgzip -d *.gz*')
    
    print()
    print('\033[1;45mYour path to reference chromosome ' + number + ' is:\033[0;0;0m')
    os.system('realpath *.fa')
        
start()
