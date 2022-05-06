#!/usr/bin/env python3
import os
import shutil

#you need to install pytul using your terminal "pip3 install python-util"
from pyutil import filereplace

#THIS PROGRAM IS FOR TRIMMED FILES ONLY

print('Use the command "readlink -f file.txt" in order to get the complete path.')

#add path to where SRA file reads are located
reads = input('Copy and paste the complete path to the directory containing your SRA files: ')
filereplace('commands_srr2.txt', 'reads_path', reads)

#add the path to where TruSeq3 file is found
tru_seq = input('Copy and paste the complete path to your TruSeq3 file: ')
filereplace('commands_srr2.txt', 'truseq3_path', tru_seq)

#add the path to where bowtie files are found (must end in "bowtie/bowtie")
bowtie = input('Copy and paste the complete path to your bowtie files: ')
filereplace('commands_srr2.txt', 'bowtie2_path', bowtie)

#add the path to where reference chromosome files are found
ref_chrom = input('Copy and paste the complete path to your reference chromosome: ')
filereplace('commands_srr2.txt', 'ref_chrom', ref_chrom)

#this asks user to type in 10 accession numbers
srr_one = input('Please paste the 1st SRA accession number: ')
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
filereplace('commands_srr2.txt',"SRR_one", srr_one)
filereplace('commands_srr2.txt',"SRR_two", srr_two)
filereplace('commands_srr2.txt',"SRR_thr", srr_thr)
filereplace('commands_srr2.txt',"SRR_fou", srr_fou)
filereplace('commands_srr2.txt',"SRR_fiv", srr_fiv)
filereplace('commands_srr2.txt',"SRR_six", srr_six)
filereplace('commands_srr2.txt',"SRR_sev", srr_sev)
filereplace('commands_srr2.txt',"SRR_eig", srr_eig)
filereplace('commands_srr2.txt',"SRR_nin", srr_nin)
filereplace('commands_srr2.txt',"SRR_ten", srr_ten)

#this will run the commands on the commands_srr2.txt file
os.system('cat commands_srr2.txt | bash')

print('*** The previous input sequences have been analyzed.')

# getting the current working directory
src_dir = os.getcwd()

# printing current directory
print('*** This is your current working directory:' + src_dir)

# printing the list of new files
print('*** These are the files in the present directory: ') 
print(os.listdir())

# copying the files
shutil.copyfile('commands_srr_template_gen2.txt', 'commands_srr2.txt') #copy src to destin

print('*** fastq_dump_tools is ready to run.')

print('*** Would you like to run another analysis? ')

while True:
    a = input('Enter yes/no to continue: ')
    if a=="yes":
        os.system('python3 fastq_dump_tools.py')
        continue
    elif a=="no":
        print('Analysis terminated.')
        break
    else:
        print('Enter either yes/no: ')



