import pandas as pd
from decimal import *
import numpy as np
import os
import shutil

#The program should be in the same directory with all combined vcf reports. It transfers them to CSV format and eliminates the SNPs with <70% in the cohort.

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# get the current working directory
cwd = os.getcwd()

# list all directories in the current working directory
dir_list = os.listdir(cwd)

# loop through the directories and find the first one that ends with "_vcf"
vcf_dir = None
for d in dir_list:
    if os.path.isdir(d) and d.endswith("_vcf"):
        # set vcf_dir to the first directory that ends with "_vcf"
        vcf_dir = d
        break

if vcf_dir is None:
    print(f"{RED}No directory ending with '_vcf' found in current directory.{RESET}")
else:
    # change the current working directory to the first directory that ends with "_vcf"
    os.chdir(vcf_dir)
    print(f"{MAGENTA}\nChanged current directory to:{RESET} ", os.getcwd())

print('\n')

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

print(f'{MAGENTA}Here is/are the vcf file(s) in your current directory:{RESET}')
num = 0
transfer = []
for i in os.listdir():
    if '.vcf' in i:
        print(i)
        num = num + 1
        transfer.append(i)
print(f'{MAGENTA}\nThere is/are{RESET} {num} {MAGENTA}vcf file(s) in your current directory.{RESET}')

vcfrun = input(f'{MAGENTA}Enter how many "comb.vcf" files you wish to convert to csv (type "ALL" to transfer all): {RESET}')

if vcfrun == 'ALL':
    for file in transfer:
        snp_above_70(file)
        print(f'{MAGENTA} {file} converted to csv format {RESET}')


if vcfrun != 'ALL':
    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
    vcfrun = int(vcfrun)
    vcfrun+=1
    placement = [ordinal(n) for n in range(1, vcfrun)]

    vcfs = []

    for time in range(int(vcfrun)-1):
        filename = input("Enter the " + placement[time] + " vcf file you wish to convert (type q to quit): ")
        if filename == 'q':
            break
        else: vcfs.append(filename)

    for name in vcfs:   
        snp_above_70(name)
        print(f'{MAGENTA}{name} converted to csv format {RESET}')

print('\n')

homedir = os.getcwd()

if 'csv_files' not in os.listdir():
        os.makedirs(os.path.join(homedir, 'csv_files'))
        print(f'{MAGENTA}Creating "csv" directory...{RESET}')

for file in os.listdir():
    if '.csv' in file:
        dst = os.path.join(homedir, 'csv_files')
        shutil.move(os.path.join(homedir, file),
                    os.path.join(dst, file))
print(f'{MAGENTA}All vcf files have been moved into the "csv" directory. {RESET}')