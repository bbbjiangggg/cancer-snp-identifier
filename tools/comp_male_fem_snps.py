import pandas as pd

# Get filenames from user
file1 = input('Enter the name of the first CSV file: ')
file2 = input('Enter the name of the second CSV file: ')

# Read CSV files
try:
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    
    # Check if "Pos" column exists in both dataframes
    if 'Pos' in df1.columns and 'Pos' in df2.columns:
        # Merge the dataframes on the "Pos" column, including all rows and adding suffixes
        merged_df = pd.merge(df1, df2, on='Pos', suffixes=('_file1', '_file2'), how='outer', indicator=True)
        
        # Find the matching rows
        matching_rows = merged_df['_merge'] == 'both'

        # Find the rows unique to each file
        unique_file1_rows = merged_df['_merge'] == 'left_only'
        unique_file2_rows = merged_df['_merge'] == 'right_only'

        # Create new dataframes for matching rows and unique rows
        matching_df = merged_df[matching_rows].drop(columns=['_merge'])
        unique_file1_df = merged_df[unique_file1_rows].drop(columns=['_merge'])
        unique_file2_df = merged_df[unique_file2_rows].drop(columns=['_merge'])

        # Write to new CSV files
        matching_df.to_csv('matching_lines.csv', index=False)
        unique_file1_df.to_csv('unique_file1.csv', index=False)
        unique_file2_df.to_csv('unique_file2.csv', index=False)

        print("Matching lines written to matching_lines.csv")
        print("Unique lines from the first file written to unique_file1.csv")
        print("Unique lines from the second file written to unique_file2.csv")
    else:
        print('"Pos" column not found in one or both of the files.')
except FileNotFoundError:
    print('One or both of the files were not found.')
except Exception as e:
    print(f"An error occurred: {str(e)}")

