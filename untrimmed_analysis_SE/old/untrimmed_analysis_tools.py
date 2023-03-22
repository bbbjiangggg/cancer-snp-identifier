#!/usr/bin/env python3
import os
import shutil


#THIS PROGRAM IS FOR WHOLE ANALYSIS

#you need to install pytul using your terminal "pip3 install python-util"
from pyutil import filereplace

#must have sendemail installed on terminal
#for ubuntu use: $ sudo apt-get install libio-socket-ssl-perl libnet-ssleay-perl sendemail
# for Mac, use: brew install sendemail

# getting the current working directory
src_dir = os.getcwd()
print('\n')
# printing current directory
print('\033[1;45m This is your current working directory:\033[0m' + src_dir + '\n')
print('Use the command "realpath filename.txt" to get the complete path. \n')

#add the email to be notified when the process is done
user = input('\033[1;45m 1) \033[0m Enter the email address to be notified once the analysis is complete: ')


#add the job title
job = input('\033[1;45m 2) \033[0m Enter a job name: ')


#add the path to where trimmomatic-0.39.jar is found
trim = input('\033[1;45m 3) \033[0m Copy and paste the complete path to your trimmomatic-0.39.jar file: ')
filereplace('untrimmed_bash_srr.txt',"trim_path", trim)

#add the path to where TruSeq3 file is found
tru_seq = input('\033[1;45m 4) \033[0m Copy and paste the complete path to your TruSeq3 file: ')
filereplace('untrimmed_bash_srr.txt', 'truseq3_path', tru_seq)

#add the path to where bowtie files are found (must end in "bowtie")
bowtie = input('\033[1;45m 5) \033[0m Copy and paste the complete path to your Bowtie files: ')
filereplace('untrimmed_bash_srr.txt', 'bowtie2_path', bowtie)

#add the path to where reference chromosome is found
ref_chrom = input('\033[1;45m 6) \033[0m Copy and paste the complete path to your BWA reference chromosome: ')
filereplace('untrimmed_bash_srr.txt', 'ref_chrom', ref_chrom)

#make a copy of untrimmed bash srr
os.system('cp untrimmed_bash_srr.txt copy_untrimmed_bash_srr.txt')


#this asks user to type in the path to the accession list
accession = input('\033[1;45m 7) \033[0m Copy and paste the name of the accesion list file (make sure it is in the same directory): ')
with open(accession,'r') as file:
    x= file.read()
    line = x.split('\n')
    line.remove('')
    print(f'There are {len(line)} unanalyzed sequences')
    file.close()

#this creates a copy of the accession list
os.system('cp ' + accession + ' copy_' + accession)
    
#this asks user to type in accession numbers
number = int(input('How many SRA sequences do you wish to analyze: '))

#this will store placement numbers into a list
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
number = number + 1
placement = [ordinal(n) for n in range(1, number)]

#this will set different variables for different srr sequences
srr_list = line[0:number-1]

 
#these commands will replace each SRR number on .txt file with 
#each of the accession numbers entered by user
for index in range(len(srr_list)):
    filereplace('untrimmed_bash_srr.txt',"number", placement[index])
    filereplace('untrimmed_bash_srr.txt',"SRR_one", srr_list[index])

    #run the commands on the untrimmed_bash_srr.txt file
    os.system('cat untrimmed_bash_srr.txt | bash')

    #replace the changed names back to orginal
    filereplace('untrimmed_bash_srr.txt', placement[index], "number")
    filereplace('untrimmed_bash_srr.txt', srr_list[index], "SRR_one")


#reset the bowtie2_path, refchrome, trimpath, and truseq path
filereplace('untrimmed_bash_srr.txt', bowtie, 'bowtie2_path')
filereplace('untrimmed_bash_srr.txt', ref_chrom, 'ref_chrom')
filereplace('untrimmed_bash_srr.txt', trim, "trim_path")
filereplace('untrimmed_bash_srr.txt', tru_seq, 'truseq3_path')

#this will remove the analyzed sequences from the accession list
with open(accession,'w') as file:
    for content in line:
        if content not in srr_list:
            file.write(content + '\n')

os.system('rm copy_untrimmed_bash_srr.txt')

#run the commands on the sendemail.txt file
print('Sending email to ' + user + ' ....')
os.system('sendemail -f sudoroot1775@outlook.com -t ' + user + ' -u ' + job + '_name_Analysis Done -m "Ready to receive information for the next analysis." -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')


print('\033[1;45m untrimmed_analysis_tools is ready to run.\033[0;0;0m \n')

print('\033[1;45m Would you like to run another analysis? \033[0;0;0m')

while True:
    a = input('Enter yes/no to continue: ')
    if a=="yes":
        print('\033[1;45m 1) \033[0m Email address used: ', user)
        print('\033[1;45m 2) \033[0m trimmomatic-0.39.jar file path: ', trim)
        print('\033[1;45m 3) \033[0m TruSeq3 file path: ', tru_seq)
        print('\033[1;45m 4) \033[0m Bowtie files path: ', bowtie)
        print('\033[1;45m 5) \033[0m BWA reference chromosome path: ', ref_chrom)
        print('\033[1;45m 6) \033[0m Accession list file name: ', accession)
        os.system('python3 untrimmed_analysis_tools.py')
        continue
    elif a=="no":
        print('\033[1;45m Analysis terminated. Goodbye. \033[0;0;0m')
        break
    else:
        print('\033[1;45m Enter either yes/no:\033[0;0;0m ')
