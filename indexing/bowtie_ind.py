#!/usr/bin/env python3

import os

#bowtie is used to index the reference chromosome

def start():
    #make necessary directory
    number = input('\033[1;45m Please enter a valid chromosome number: \033[0;0;0m')
    print('\n')

    isvalid = number.startswith('0')
    greater = int(number) < 23

    while ((isvalid == False) and (greater == True)):
        os.system('mkdir ' + number + '_bowtie_ind')
        os.chdir(number + '_bowtie_ind')
        
        #download chromosome sequence
        print('\033[1;45m Downloading Chromosome '+ number + '\033[0;0;0m' )
        os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.'+ number +'.fa.gz')
      
        #indexing chromosome file
        print('\033[1;45m Indexing with Bowtie... \033[0;0;0m' )
        os.system('bowtie2-build *.gz* bowtie')
        
        print('\033[1;45m Indexing with Bowtie complete \033[0;0;0m' )
        print('\n')
        print('\033[1;45m Your path to the Bowtie index of chromsome ' + number + ' is: \033[0;0;0m' )
        os.system('realpath bowtie')
        exit()
    else:
        print('\033[0;101m NOTE: Do not include a "0" in front of the number. \033[0;0;0m')
        print('\033[0;101m The number must be less than 23. \033[0;0;0m')
        
        print('\n')
        start()
     

start()
    

            