import os
from os import listdir
from os.path import join, isfile
import shutil
from shutil import copyfile


print('Use the command "readlink -f name_of_file/dir" in order to get the complete path.')
print("\033[1;35;40m Bright Magenta \033[0m 1;35;40m")

directory = input('\033[1;35;40m Enter the name of the directory you wish to create and copy all vcf files to: \033[0m 0;0;0m')

homedir = input('Enter the path where the SRR directories are currently stored: ')

isecdir = os.path.join(homedir, directory)

os.mkdir(isecdir)
print('\033[1;35;40m New directory has been created: \033[0m 0;0;0m' + directory)
print('\033[1;35;40m Copying all vcf files to: ' + directory + ' \033[0m 0;0;0m')

for item in listdir(homedir):
    if 'SRR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                copyfile(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('\033[1;35;40m All vcf files have been copied to: ' + directory + ' \033[0m 0;0;0m')

