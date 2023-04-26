import csv

# Function to count number of 1s in a row of column 5
def count_ones(row):
    return row.count('1')

# Prompt user for input file path
input_file_path = input("Enter the input VCF file path: ")

# Prompt user for output file path
output_file_path = input("Enter the output CSV file path: ")

# Open input file and read contents
with open(input_file_path, 'r') as f:
    lines = f.readlines()

# Loop through each line in file, skipping header lines
updated_lines = []
for line in lines:
    if line.startswith('#'):
        continue
    
    # Split line into fields and count number of 1s in column 5
    fields = line.split('\t')
    col5 = fields[4]
    num_ones = count_ones(col5)

    # Add number of ones as last column and append updated line
    fields.append(str(num_ones))
    updated_lines.append(fields)

# Write updated lines to output file in CSV format
with open(output_file_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(updated_lines)

# Print confirmation message
print("VCF file converted to CSV format and written to output file:", output_file_path)
