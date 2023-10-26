# Prompt the user to enter the chromosomes
chromosomes_input = input("Enter the chromosome numbers separated by comma (e.g., 1,2,X,Y) to be analyzed: ")

# Split the input string by comma
chromosomes_list = chromosomes_input.split(',')

# Remove leading and trailing whitespace from each chromosome name
chromosomes_list = [chromosome.strip() for chromosome in chromosomes_list]

# Print the list of chromosomes
print("List of chromosomes to be analyzed:", chromosomes_list)
