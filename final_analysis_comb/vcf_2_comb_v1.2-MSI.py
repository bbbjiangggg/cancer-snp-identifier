def count_ones_in_vcf(file_path):
    # Open the input VCF file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    output_lines = []

    # Process each line in the VCF file
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            # Skip header lines starting with '#'
            output_lines.append(line)
        else:
            fields = line.split('\t')
            info = fields[7]  # Fifth column (index 4) in VCF format

            # Count the number of ones in the info field
            ones_count = info.count('1')

            # Append the count as a new column
            output_line = f'{line}\t{ones_count}'
            output_lines.append(output_line)

    # Write the output to a new VCF file
    output_file_path = 'output.vcf'
    with open(output_file_path, 'w') as output_file:
        output_file.write('\n'.join(output_lines))

    print(f"Counting completed. The results are saved in '{output_file_path}'.")


# Prompt the user to enter the path of the merged VCF file
merged_vcf_file = input("Enter the path of the merged VCF file: ")

# Call the function to process the VCF file
count_ones_in_vcf(merged_vcf_file)
