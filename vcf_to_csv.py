import pandas as pd
from decimal import *
import numpy as np

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
        #writer = pd.ExcelWriter(filename + '.xlsx', engine='xlsxwriter')
        #sort.to_excel(writer, sheet_name='welcome', index=False)
        #writer.save()
        print(sort)

vcfrun = input('Enter how many vcf file do you wish to transfer to csv file: ')

for time in range(int(vcfrun)):
    filename = input("Enter the " + str(time+1) + " vcf file you wish to analyze (type q to quit): ")
    if filename == 'q':
        break
    snp_above_70(filename)


