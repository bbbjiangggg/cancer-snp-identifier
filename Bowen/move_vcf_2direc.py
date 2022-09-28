import os
from os import listdir
from os.path import join, isfile
import shutil

#Copying all vcf reports to one directory

print('*** Use the command "readlink -f name_of_file/dir" in order to get the complete path.')

directory = input('>>> Enter the name of the directory you wish to create and copy all vcf files to (e.g. isec_prca_ch1):\n')

homedir = input('>>> Enter the path where the SRR directories are currently stored: ')

isecdir = os.path.join(homedir, directory)

os.mkdir(isecdir)
print('\033[1;45m New directory has been created: \033[0m' + directory)
print('\033[1;45m Copying all vcf files to: ' + directory + ' \033[0m')

for item in listdir(homedir):
    if 'SRR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                shutil.copyfile(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('\033[1;35;40m All vcf files have been copied to: ' + directory + ' \033[0m 0;0;0m')

copyName = 'copy_'+directory
copydir = os.path.join(homedir, copyName)
shutil.copytree(isecdir, copydir)

print('\033[1;35;40m A copy of ' + directory + ' has been made \033[0m')

#Combining the vcf reports

name = input('\033[1;42m Enter the chromosome and cancer name for you analysis with a dash in between (e.g. ch10_prca): \033[0m')

os.chdir(os.path.join(homedir, isecdir))

print('\033[1;45m Compressing all vcf files to gzip... \033[0;0;0m')
os.system('ls *.vcf | xargs -n1 -P0 bgzip')
print('\n')

os.system('ls -1 | paste -sd " " ->> isec_commands.txt')

print('\033[1;45m Indexing all the gzip files \033[0;0;0m')
os.system('for f in ./*.vcf.gz; do tabix -p vcf -f $f;done')
print('\n')

print('\033[1;45m Modifying the text commands \033[0;0;0m')
with open('isec_commands.txt', 'r+') as f:
    text = f.read()
    text = text.replace('isec_commands.txt', '| bgzip -c >bgzip -c > ' + name + '_comb.vcf.gz')
    f.seek(0, 0)
    f.write('bcftools isec -n +2' + ' ' + text) 
    f.close()
print('\n')

print('\033[1;45m Running the text commands \033[0;0;0m')
os.system('cat isec_commands.txt | bash')
print('\n')

print('\033[1;45m Unzipping the combined vcf report \033[0;0;0m')
os.system('gunzip ' + name + '_comb.vcf.gz')
print('\n')

print('\033[1;45m Final vcf report created \033[0;0;0m')

