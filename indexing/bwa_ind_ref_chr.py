#!/usr/bin/env python3

import os

#BWA is used to index the reference chromosome

def start():
    #make necessary directory
    number = input('\033[1;45m Please enter a valid chromosome number: \033[0;0;0m')
    print('\n')

    isvalid = number.startswith('0')
    greater = int(number) < 23

    while ((isvalid == False) and (greater == True)):
        os.system('mkdir bwa_ind_'+ number)
        os.system('cd bwa_ind_'+ number)
        
        print('\033[1;45m Downloading Chromosome '+ number + '\033[0;0;0m' )
        os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.'+ number +'.fa.gz')

        #download chromosome sequence
        os.system('mv *.gz* bwa_ind_'+ number)

        #indexing chromosome file
        print('\033[1;45m Indexing with BWA... \033[0;0;0m' )
        os.chdir('bwa_ind_' + number)
        os.system('bwa index *.gz*')
        os.system('bgzip -d *.gz*')
        print('\033[1;45m Indexing with BWA complete \033[0;0;0m' )
        print('\n')
        print('\033[1;45m Your path to reference chromosome ' + number + ' is: \033[0;0;0m' )
        os.system('realpath *.fa')
        exit()
    else:
        print('\033[0;101m NOTE: Do not include a "0" in front of the number. \033[0;0;0m')
        print('\033[0;101m The number must be less than 23. \033[0;0;0m')
        
        print('\n')
        start()
     

start()
    

            