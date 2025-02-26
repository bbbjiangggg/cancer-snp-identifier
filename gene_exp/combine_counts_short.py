import os
import pandas as pd

# Find all CSV files in the current directory
csv_files = [file for file in os.listdir() if file.endswith('.csv')]

# Combine the CSV files into a single DataFrame and average the counts
combined_df = pd.concat([pd.read_csv(file, index_col=0) for file in csv_files], axis=1)

# Calculate the average counts per gene
combined_df['average_count'] = combined_df.mean(axis=1)

# Drop rows with specified index values
combined_df = combined_df.drop(['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique'], errors='ignore')

# Sort by average count in descending order
sorted_df = combined_df[['average_count']].sort_values(by='average_count', ascending=False)

# Write only the Gene ID and average count to the output file
sorted_df.to_csv('count_comb_output.csv', header=['average_count'])

print("Successfully created 'count_comb_output_short.csv' with Gene ID and Average Counts!")
