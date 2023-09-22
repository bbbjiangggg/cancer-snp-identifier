import os
import re
from os import listdir
from os.path import join, isfile
import shutil
from shutil import copyfile


print('Use the command "readlink -f name_of_file/dir" in order to get the complete path.')


directory = input('\033[1;45m Enter the name of the directory you wish to create and copy all vcf files to: \033[0;0;0m ')

print('\n')
homedir = input('\033[1;45m Enter the path where the SRR directories are currently stored: \033[0;0;0m ')


isecdir = os.path.join(homedir, directory)

os.mkdir(isecdir)
print('\n')
print('\033[1;45m New directory has been created: \033[0;0;0m ' + directory)

print('\n')
print('\033[1;45m Copying all vcf files to: \033[0;0;0m ' + directory)


for item in listdir(homedir):
    if 'SRR' or 'ERR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                copyfile(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass
print('\n')
print('\033[1;45m All vcf files have been copied to: \033[0;0;0m ' + directory)

copyName = 'copy_'+directory
copydir = os.path.join(homedir, copyName)
shutil.copytree(isecdir, copydir)

print('\n')
print('\033[1;45m A copy of the directory has been created: \033[0;0;0m ' + directory)

print('\n')

#combining all vcf files into one

combo = input('\033[1;45m Enter the chromosome and cancer type (e.g. ch2_pnca): \033[0;0;0m ')
os.chdir(os.path.join(homedir, isecdir))

print('\033[1;45m Compressing all vcf files to gzip... \033[0;0;0m')
os.system('ls *.vcf | xargs -n1 -P0 bgzip')
print('\n')

#this command will Write all file names into a txt file, same line, 
#one space, named "isec_tools_commands.txt"
os.system('ls -1 | paste -sd " " ->> isec_tools_commands.txt')
print('\n')

print('\033[1;45m Indexing vcf files using tabix... \033[0;0;0m')

#tabix all vcfgz files
os.system('ls *.gz | xargs -n1 -P0 tabix')
print('\n')

print('\033[1;45m Modifying the text commands \033[0;0;0m')
with open('isec_commands.txt', 'r') as file:
    data = file.read()
    data = data.replace('isec_commands.txt', ' | bgzip -c >bgzip -c > ' + combo + '_comb.vcf.gz')
    file.seek(0, 0)
    file.write('bcftools isec -n=2 -p ' + combo + '_comb ' + data)
    file.close()
print('\n')

print('\033[1;45m Running the isec tools commands... \033[0;0;0m')
os.system('bash isec_commands.txt')
print('\n')

print('\033[1;45m Indexing the combined vcf file... \033[0;0;0m')
os.system('gunzip ' + combo + '_comb.vcf.gz')
print('\n')

print('\033[1;45m Final vcf report created \033[0;0;0m')