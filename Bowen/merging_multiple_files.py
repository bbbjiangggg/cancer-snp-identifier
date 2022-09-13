import pandas as pd
from functools import reduce
import os

#The program combines CSV files for different cancers and finds the common SNPs

dataCount = 0
data_frames = []
csvs = []

print('**********************************************************************')
print('*****  Here are the csv files in your current direcotry: ')
num = 0
for i in os.listdir():
    if '.csv' in i:
        csvs.append(i)
        print(i)
        num += 1
print(f'***** You have {num} csv files in your current directory')

repeats = input('>>> How many csv files do you wish to combine (Enter a number or "ALL" to merge all): ')
chromosome = input('>>> What is the chromosome that you are analysing (Enter like ch1 or ch10): ')


if repeats == 'ALL':
    for repeat in range(num):
        globals()['df%s'%(repeat+1)] = pd.read_csv(csvs[repeat], index_col = [0])
        globals()['df%s'%(repeat+1)].pop('Individual')

        file = csvs[repeat]
        
        filename = file.replace('.csv', '')
        filename = filename.replace(chromosome, '')
        filename = filename.replace('comb', '')
        filename = filename.replace('_', '')

        globals()['df%s'%(repeat+1)].rename(columns={"Chrom": ("Chrom_" + str(filename)), "Ref": ("Ref_" + str(filename)),
                                                 "Mut": ("Mut_" + str(filename)), "Count": ("Count_" + str(filename))}, inplace = True)

        dataCount = dataCount + 1
        
if repeats != 'ALL':
    for repeat in range(int(repeats)):
        file = input('> Enter the '+ str(repeat+1) + ' csv file: ')

        filename = file.replace('.csv', '')
        filename = filename.replace(chromosome, '')
        filename = filename.replace('comb', '')
        filename = filename.replace('_', '')
        
        globals()['df%s'%(repeat+1)] = pd.read_csv(file, index_col = [0])
        globals()['df%s'%(repeat+1)].pop('Individual')
        
        globals()['df%s'%(repeat+1)].rename(columns={"Chrom": ("Chrom_" + str(filename)), "Ref": ("Ref_" + str(filename)),
                                                     "Mut": ("Mut_" + str(filename)), "Count": ("Count_" + str(filename))}, inplace = True)
        
        dataCount = dataCount + 1

for num in range(dataCount):
    data_frames.append(globals()['df%s'%(num+1)])

df_merged = reduce(lambda left,right: pd.merge(left, right, on=['Pos'], how='outer'), data_frames)

nan_value = float("NaN")
df_merged.replace("", nan_value, inplace=True)

for column in df_merged.columns.tolist():
    df_merged.dropna(subset = [column], inplace=True)

wish = input('>>> What do you want to name your merged file (add .csv in the end )? ')
df_merged.to_csv(wish)

print('\033[1;45m Merging Completed \033[0m')
