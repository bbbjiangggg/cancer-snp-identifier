import pandas as pd

# Ask the user to input the CSV file
csv_file = input("Enter the path to the CSV file: ")

# Read the CSV file into a DataFrame
data = pd.read_csv(csv_file, delimiter=',')

# Print the column names
print("Column Names:")
print(data.columns)

# Ask the user to input the column name for sorting
sort_column = input("Enter the column name for sorting: ")

# Sort the DataFrame by the specified column in descending order
sorted_data = data.sort_values(by=sort_column, ascending=False)

# Get the file name without the extension
file_name = csv_file.split(".")[0]

# Create the output file name
output_file = "sorted_" + file_name + ".csv"

# Write the sorted data to a new CSV file
sorted_data.to_csv(output_file, index=False)

# Print a success message
print("Sorted data has been written to", output_file)

