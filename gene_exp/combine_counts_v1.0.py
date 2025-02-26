import os
import pandas as pd
from termcolor import colored

# Find all CSV files in the current directory
print(colored("🔍 Searching for CSV files...", "cyan"))
csv_files = [file for file in os.listdir() if file.endswith('.csv')]

if not csv_files:
    print(colored("⚠️ No CSV files found in the current directory!", "red"))
    exit()

print(colored(f"📂 Found {len(csv_files)} CSV files. Processing...", "cyan"))

# Combine the CSV files into a single DataFrame and average the counts
print(colored("📊 Reading and merging CSV files...", "yellow"))
combined_df = pd.concat([pd.read_csv(file, index_col=0) for file in csv_files], axis=1)

# Calculate the average counts per gene
print(colored("🧬 Calculating average counts per gene...", "yellow"))
combined_df['average_count'] = combined_df.mean(axis=1)

# Drop unwanted rows
print(colored("🗑️ Removing non-gene entries...", "yellow"))
combined_df = combined_df.drop(['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique'], errors='ignore')

# Sort by average count in descending order
print(colored("🔢 Sorting genes by average count...", "yellow"))
sorted_df = combined_df[['average_count']].sort_values(by='average_count', ascending=False)

# Write only the Gene ID and average count to the output file
output_file = "count_comb_output.csv"
sorted_df.to_csv(output_file, header=['average_count'])

print(colored(f"✅ Successfully created '{output_file}' with Gene ID and Average Counts!", "green"))
