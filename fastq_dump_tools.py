#!/usr/bin/env python3
import os

#you need to install pytul using your terminal "pip3 install python-util"
from pyutil import filereplace

#add the path to where trimmomatic-0.39.jar is found
trim = input('Copy and paste the complete path to your trimmomatic-0.39.jar file: ')
filereplace('commands_srr.txt',"trim_path", trim)

#add the path to where TruSeq3 file is found
tru_seq = input('Copy and paste the complete path to your TruSeq3 file: ')
filereplace('commands_srr.txt', 'truseq3', tru_seq)

#add the path to where bowtie files are found (must end in "bowtie/bowtie")
bowtie = input('Copy and paste the complete path to your bowtie files: ')
filereplace('commands_srr.txt', 'bowtie2_path', bowtie)

ref_chrom = input('Copy and paste the complete path to your reference chromosome: ')
filereplace('commands_srr.txt', '/mnt/d/sra/brca/bgz/Homo_sapiens.GRCh38.dna.chromosome.10', ref_chrom)

#this asks user to type in 10 accession numbers
'''srr_one = input('Please paste the 1st SRA accession number: ')
srr_two = input('Please paste the 2nd SRA accession number: ')
srr_thr = input('Please paste the 3rd SRA accession number: ')
srr_fou = input('Please paste the 4th SRA accession number: ')
srr_fiv = input('Please paste the 5th SRA accession number: ')
srr_six = input('Please paste the 6th SRA accession number: ')
srr_sev = input('Please paste the 7th SRA accession number: ')
srr_eig = input('Please paste the 8th SRA accession number: ')
srr_nin = input('Please paste the 9th SRA accession number: ')
srr_ten = input('Please paste the 10th SRA accession number: ')



#these commands will replace each SRR number on .txt file with 
#each of the ten accession numbers entered by user
filereplace('commands_srr.txt',"SRR_one", srr_one)
filereplace('commands_srr.txt',"SRR_two", srr_two)
filereplace('commands_srr.txt',"SRR_thr", srr_thr)
filereplace('commands_srr.txt',"SRR_fou", srr_fou)
filereplace('commands_srr.txt',"SRR_fiv", srr_fiv)
filereplace('commands_srr.txt',"SRR_six", srr_six)
filereplace('commands_srr.txt',"SRR_sev", srr_sev)
filereplace('commands_srr.txt',"SRR_eig", srr_eig)
filereplace('commands_srr.txt',"SRR_nin", srr_nin)
filereplace('commands_srr.txt',"SRR_ten", srr_ten)

#this will run the commands on the commands_srr.txt file
os.system('cat commands_srr.txt | bash')'''
