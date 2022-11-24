import os
from os import listdir
from os.path import join, isfile
import shutil
import re
import pandas as pd
from decimal import *
import numpy as np

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

#The program should be in the same directory with all combined vcf reports. It transfers them to CSV format and eliminates the SNPs with <70% in the cohort.

#os.chdir(directory + '/isec_vcfgz_files/')

def snp_above_70(filename):
    getcontext().prec = 2
    count = 0
    individual_count = 0
    allcols = []

    #check the count percentage
    with open(filename) as file:
        for line in file:
            cols = line.split('\t')
            for letter in cols[4]:
                individual_count = individual_count + 1
                if letter == "1":
                    count += 1
            if Decimal(count)/Decimal(individual_count) >= 0.70:            
                cols.append(str(count))
            allcols.append(cols)
            count = 0
            individual_count = 0


        #delete rows with NaN and sort the other counts
        df = pd.DataFrame(allcols)
        df.columns = ['Chrom', 'Pos', 'Ref', 'Mut', 'Individual', 'Count']
        df.Count = pd.to_numeric(df.Count, errors='coerce')
        df['Count'].replace('', np.nan, inplace=True)
        df.dropna(subset=['Count'], inplace=True)
        sort = df.sort_values(by = 'Count', ascending=False)
        filename = filename.replace('.vcf', '')
        sort.to_csv(filename+'.csv')
        print(sort)

print('\033[1;45mMake sure your combined vcf file(s) is/are in the same directory with this program.\033[0;0;0m')
print('\033[1;45mHere is/are the vcf file(s) in your current directory:\033[0;0;0m')
num = 0
transfer = []
for i in os.listdir():
    if '.vcf' in i:
        print(i)
        num = num + 1
        transfer.append(i)
print(f'\033[1;45mThere is/are {num} vcf file(s) in your current directory.\033[0;0;0m')

vcfrun = input('\033[1;45m Enter how many vcf files you wish to transfer to the csv directory (type "ALL" to transfer all): \033[0m')

if vcfrun == 'ALL':
    for file in transfer:
        snp_above_70(file)
        print(f'\033[1;45m {file} changed to csv format \033[0m')


if vcfrun != 'ALL':
    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
    vcfrun = int(vcfrun)
    vcfrun+=1
    placement = [ordinal(n) for n in range(1, vcfrun)]

    vcfs = []

    for time in range(int(vcfrun)-1):
        filename = input("Enter the " + placement[time] + " vcf file you wish to analyze (type q to quit): ")
        if filename == 'q':
            break
        else: vcfs.append(filename)

    for name in vcfs:   
        snp_above_70(name)
        print(f'\033[1;45m {name} changed to csv format \033[0m')

homedir = os.getcwd()

if 'CSV_files' not in os.listdir():
        os.makedirs(os.path.join(homedir, 'CSV_files'))

for file in os.listdir():
    if '.csv' in file:
        dst = os.path.join(homedir, 'CSV_files')
        shutil.move(os.path.join(homedir, file),
                    os.path.join(dst, file))

print('\033[1;45m All vcf files have been moved into the CSV directory. \033[0m')