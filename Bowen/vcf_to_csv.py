import pandas as pd
from decimal import *
import numpy as np
import os
import shutil

#The program should be in the same directory with all combined vcf reports. It transfers them to CSV format and eliminates the SNPs with <70% in the cohort.

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

print('***** Make sure your vcf files are in the same directory with the code *****')
print('***** Here are the vcf files in your current directory: ')
num = 0
transfer = []
for i in os.listdir():
    if '.vcf' in i:
        print(i)
        num = num + 1
        transfer.append(i)
print(f'***** There are {num} vcf files in your current directory *****')

vcfrun = input('\033[1;45m Enter how many vcf file do you wish to transfer to csv file (type "ALL" to transfer all): \033[0m')

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

print('\033[1;45m All vcf files are moved into the CSV folder \033[0m')

