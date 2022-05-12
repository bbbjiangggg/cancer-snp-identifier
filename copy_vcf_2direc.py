import os
from os import listdir
from os.path import join, isfile
import shutil
from shutil import copyfile


print('*** Use the command "readlink -f name_of_file/dir" in order to get the complete path. ***')

directory = input('Enter the name of the directory you wish to create and copy all vcf files to: ')

homedir = input('Enter the path to the directory where the SRR directories are currently stored: ')

isecdir = os.path.join(homedir, directory)

os.mkdir(isecdir)
print('*** New directory has been created: ' + directory + ' ***')
print('*** Copying all vcf files to: ' + directory + ' ***')

for item in listdir(homedir):
    if 'SRR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                copyfile(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('*** All vcf files have been copied to: ' + directory + ' ***')