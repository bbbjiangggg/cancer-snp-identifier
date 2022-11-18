#!/usr/bin/env python3

import os

#bowtie is used to index the reference chromosome

def start():
    #make necessary directory
    number = input('\033[1;45m Please enter chromosome X or Y: \033[0;0;0m')
    print('\n')

    isvalid = number.startswith('X')
    isvalid2 = number.startswith('Y')
    

    while ((isvalid == True) or (isvalid2 == True)):
        if number.startswith('X'):
            os.system('mkdir bowtie_ind_'+ number)
            os.chdir('bowtie_ind_' + number)
            print('\033[1;45m Downloading chromosome ' + number + ' sequence \033[0;0;0m' )
            os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.'+ number +'.fa.gz')
            
            #indexing chromosome file
            print('\033[1;45m Indexing with Bowtie... \033[0;0;0m' )
            os.system('bowtie2-build *.gz* bowtie')
            print('\033[1;45m Indexing of chromosome ' + number + ' with bowtie complete \033[0;0;0m' )
            print('\n')
            print('\033[1;45m Your path to reference chromosome ' + number + ' is: \033[0;0;0m' )
            os.system('realpath bowtie')
            break
        
        elif number.startswith('Y'):
            os.system('mkdir bowtie_ind_'+ number)
            os.chdir('bowtie_ind_' + number)
            print('\033[1;45m Downloading chromosome ' + number + ' sequence \033[0;0;0m' )
            os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.'+ number +'.fa.gz')
            
            #indexing chromosome file
            print('\033[1;45m Indexing with Bowtie... \033[0;0;0m' )
            os.system('bowtie2-build *.gz* bowtie')
            print('\033[1;45m Indexing of chromosome ' + number + ' with bowtie complete \033[0;0;0m' )
            print('\n')
            print('\033[1;45m Your path to reference chromosome ' + number + ' is: \033[0;0;0m' )
            os.system('realpath bowtie')
            break
        
        
    else:
        print('\033[0;101m NOTE: You must enter a capital "X" or "Y". \033[0;0;0m')
              
        print('\n')
        start()
start()  

            