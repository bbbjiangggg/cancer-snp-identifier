import os
from os import listdir
from os.path import join, isfile
import shutil
import re

#Copying all vcf reports to one directory

print('*** Use the command "realpath name_of_file/dir" to get the complete path. ***')

chr = input('\033[1;45m Enter the chrosome number of this analysis: \033[0;0;0m')
print('\n')
cancer = input('\033[1;45m Enter an abbreviation for the cancer type (e.g. pnca): \033[0;0;0m')
print('\n')

directory = 'ch' + chr + '_' + cancer + '_vcf'
print("\033[1;45m This is your directory's name: \033[0;0;0m" + directory)
print('\n')
homedir = input('\033[1;45m Enter the path where the SRR directories are currently stored: \033[0;0;0m')
print('\n')

isecdir = os.path.join(homedir,directory)

#make new directory
if directory in os.listdir():
    print('\033[1;45m Directory already exists, removing directory. \033[0;0;0m')
    os.system('rm -r ' + directory)
    os.mkdir(directory)
else:
    os.mkdir(directory)
    print('\033[1;45m New directory has been created: \033[0;0;0m ' + directory)
    print('\n')
input('\033[1;45mPress enter to continue...\033[0;0;0m')
os.chdir(isecdir)
#copy all vcf files to new directory
print('\033[1;45m Copying all vcf files to: \033[0;0;0m ' + directory + '\n')

for item in listdir(homedir):
    if 'RR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                shutil.move(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('\033[1;35;40m All vcf files have been copied to: ' + directory + ' \033[0;0;0m')
print('\n')

copyName = 'copy_'+directory
if copyName in os.listdir():
    print('\033[1;45m Directory already exists, removing directory. \033[0;0;0m')
    os.system('rm -r ' + copyName)
    os.mkdir(copyName)
else:
    #copydir = os.path.join(homedir, copyName)
    shutil.copytree(isecdir, copyName)

print('\033[1;35;40m As backup, a copy of ' + directory + ' has been created. \033[0m')
input('\033[1;35;40m Press enter to continue...\033[0m')
print('\n')

#Combining vcf reports
print('\033[1;45m Combining all vcf files in: \033[0;0;0m ' + directory + '\n')
#this command will bgzip all .vcf files, this program must be located in the same directory
#if all files are already bgzipped then it will throw back:
#"ls: cannot access '*.vcf': No such file or directory"

print('\033[1;45m Compressing all vcf files to gzip... \033[0;0;0m')


os.system('ls *.vcf | xargs -n1 -P0 bgzip')
print('\n')
print('\033[1;45m Organizing all compressed vcfgz files... \033[0;0;0m')
print('\n')

#move all .vcf files to a directory called isec_vcfgz_files
os.system('mkdir isec_vcfgz_files')
os.system('mv *.gz* isec_vcfgz_files')

#change directory to isec_vcfgz_files
os.chdir('isec_vcfgz_files')



#this command will Write all file names into a txt file, same line, 
#one space, named "isec_tools_commands.txt"
os.system('ls -1 | paste -sd " " ->> isec_tools_commands.txt')

#deleting string that results from previous command
infile = "isec_tools_commands.txt"
outfile = "isec_tools_commands2.txt"

delete_string = ["isec_tools_commands.txt"]
fin = open(infile)
fout = open(outfile, "w+")
for line in fin:
    for word in delete_string:
        line = line.replace(word, "")
    fout.write(line)
fin.close()
fout.close()



#read the existing text from file in READ mode 
with open('isec_tools_commands2.txt','r') as src:
    fline='bcftools isec -n +2 '
    #Prepending string
    oline=src.readlines()
    #prepend the string we want to, on first line 
    oline.insert(0,fline) 
#open the file in WRITE mode  
with open('isec_tools_commands2.txt','w') as src: 
   src.writelines(oline) 

vcf_file = input('\033[1;45m Enter the name you wish to give the combined vcf files (e.g. ch22_pnca_comb): \033[0m')

#open a file with access mode 'a+'web: https://stackabuse.com/file-handling-in-python/
with open('isec_tools_commands2.txt', 'a+') as file_object:
    # Append "| bgzip -c >bgzip -c > vcf_file.vcf.gz" at the end of file
    file_object.write('| bgzip -c >bgzip -c > ' + vcf_file + '.vcf.gz')
    
#remove all white spaces > 1 and save it as "isec_tools_commands2.txt"
with open('isec_tools_commands2.txt', 'r') as file_object, open ('isec_tools_commands3.txt', 'w') as file_object2:
    for line in file_object:
        file_object2.write(re.sub('\s+',' ',line))

print('\n')
print('\033[1;45m Indexing vcf files using tabix... \033[0;0;0m')


#tabix all vcfgz files
os.system('for f in ./*.vcf.gz; do tabix -p vcf -f $f;done')
print('\n')
print('\033[1;45m Creating intersections, unions and complements (isec)... \033[0;0;0m')

#run isec on isec_tools_commands2.txt files
os.system('cat isec_tools_commands3.txt | bash')

#unzip gzip file
os.system('gunzip '+ vcf_file + '.vcf.gz')

#open unzipped file
print('\n')
print('\033[1;45m Done! Your final combined vcf file is located here:\033[0;0;0m' + directory + '/isec_vcfgz_files/' + vcf_file + '.vcf ') 
print(' \033[1;45m Open the file using the command:\033[0;0;0m cat ' + directory + '/isec_vcfgz_files/' + vcf_file + '.vcf | less')


os.system('rm bgzip isec_tools_commands.txt isec_tools_commands2.txt isec_tools_commands3.txt')

