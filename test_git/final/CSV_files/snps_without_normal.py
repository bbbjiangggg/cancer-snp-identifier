import pandas as pd

#The program eliminates the normal SNPs

merged_file = input('>>> Enter the csv merged cancer file: ')

normal = input('>>> Enter the full name of the normal SNP file (make sure the normal SNP file is in the same directory): ')

ch = input('>>> What chromosome are you analyzing (e.g. ch1)?')

dfn = pd.read_csv(normal, index_col=[0])
dfm = pd.read_csv(merged_file, index_col=[0])

positions = dfn['Pos'].tolist()

dfm = dfm.loc[dfm['Pos'].isin(positions)==False]

name = ch + '_final_merge.csv'

dfm.to_csv(name, index=False)
