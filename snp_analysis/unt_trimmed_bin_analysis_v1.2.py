# Prompt the user to enter the chromosomes
chromosomes_input = input("Please enter the chromosomes to be analyzed, separated by a comma: ")

# Split the input string by comma
chromosomes_list = chromosomes_input.split(',')

# Remove leading and trailing whitespace from each chromosome name
chromosomes_list = [chromosome.strip() for chromosome in chromosomes_list]

# Print the list of chromosomes
print("List of chromosomes to be analyzed:", chromosomes_list)

# Define the base paths for BWA and Bowtie
bwa_base_path = "/usr/local/bin/bwa/"
bowtie_base_path = "/usr/local/bin/bowtie/"

# Print the paths for each chromosome
for chromosome in chromosomes_list:
    # Construct the paths
    bwa_chrom_path = f"{bwa_base_path}{chromosome}_bwa_ind/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
    bowtie_index_path = f"{bowtie_base_path}{chromosome}_bowtie_ind/bowtie"

    # Print the paths
    print(f"\nPaths for chromosome {chromosome}:")
    print("BWA Chromosome Path:", bwa_chrom_path)
    print("Bowtie Index Path:", bowtie_index_path)
