#!/usr/bin/env python3

import os
import re

#this command will bgzip all .vcf files, this program must be located in the same directory
#if all files are already bgzipped then it will throw back:
#"ls: cannot access '*.vcf': No such file or directory"

print('\033[1;45m Zipping all vcf files \033[0;0;0m')

os.system('ls *.vcf | xargs -n1 -P0 bgzip')

#move all .vcf files to a directory called isec_vcfgz_files
os.system('mkdir isec_vcfgz_files')
os.system('mv *.gz* isec_vcfgz_files')

#change directory to isec_vcfgz_files
os.chdir('isec_vcfgz_files')


#this command will Write all file names into a txt file, same line, 
#one space, named "isec_tools_commands.txt"
os.system('ls -1 | paste -sd " " ->> isec_tools_commands.txt')

#read the existing text from file in READ mode 
with open('isec_tools_commands.txt','r') as src:
    fline='bcftools isec -n +2 '
    #Prepending string
    oline=src.readlines()
    #prepend the string we want to, on first line 
    oline.insert(0,fline) 
#open the file in WRITE mode  
with open('isec_tools_commands.txt','w') as src: 
   src.writelines(oline) 

#removes the string "isec_tools_commands.txt" from the file isec_tools_commands.txt
os.system('sed -e s/isec_tools_commands.txt//g -i * isec_tools_commands.txt')


vcf_file = input('Enter the name you wish to give the combined vcf files: ')

#open a file with access mode 'a+'web: https://stackabuse.com/file-handling-in-python/
with open('isec_tools_commands.txt', 'a+') as file_object:
    # Append "| bgzip -c >bgzip -c > ch10_prca_cfDNA_comb.vcf.gz" at the end of file
    file_object.write('| bgzip -c >bgzip -c > ' + vcf_file)
    
#remove all white spaces > 1 and save it as "isec_tools_commands2.txt"
with open('isec_tools_commands.txt', 'r') as file_object, open ('isec_tools_commands2.txt', 'w') as file_object2:
    for line in file_object:
        file_object2.write(re.sub('\s+',' ',line))

print('\033[1;45m Indexing vcf files using tabix \033[0;0;0m')

#tabix all vcfgz files
os.system('for f in ./*.vcf.gz; do tabix -p vcf -f $f;done')

print('\033[1;45m Creating intersections, unions and complements (isec) \033[0;0;0m')

#run isec on isec_tools_commands2.txt files
os.system('cat isec_tools_commands2.txt | bash')

print('\033[1;45m Opening combined vcf files \033[0;0;0m')

#unzip gzip file
os.system('gunzip '+ vcf_file)

#open unzipped file
os.system('cat ' + vcf_file + ' | less')
