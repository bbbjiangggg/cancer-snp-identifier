import pandas as pd

# Prompt the user to enter the filename
filename = input("Enter the name of the CSV file: ")

# Read the Excel file into a DataFrame
df = pd.read_csv(filename)

# Get the name of the last column
last_column = df.columns[-1]

# Count the frequency of each unique entry in the last column
value_counts = df[last_column].value_counts()

# Print the result
print(f'The last column "{last_column}" contains {len(value_counts)} unique values:')
print(value_counts)

