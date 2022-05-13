import os
from os import listdir
from os.path import join, isfile
import shutil
from shutil import copyfile


print('\033[1;45m Use the command "readlink -f name_of_file/dir" in order to get the complete path.\033[0;0;0m')

print('\n')

directory = input('Enter the name of the directory you wish to create to hold all vcf files: ')

homedir = input('Enter the path to the directory where the SRR directories are currently stored: ')

isecdir = os.path.join(homedir, directory)

os.mkdir(isecdir)
print('\033[1;45m New directory has been created: ' + directory + ' \033[0;0;0m')
print('\n')
print('\033[1;45m Copying all vcf files to: ' + directory + ' \033[0;0;0m')

for item in listdir(homedir):
    if 'SRR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                copyfile(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('\033[1;45m All vcf files have been copied to: ' + directory + ' \033[0;0;0m')