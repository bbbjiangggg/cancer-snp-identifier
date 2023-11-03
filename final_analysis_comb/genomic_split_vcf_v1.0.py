import os
import glob

def is_valid_chromosome(chromosome):
    """Check if the chromosome name is valid (1-22, X, Y)."""
    return chromosome.isdigit() and 1 <= int(chromosome) <= 22 or chromosome in ['X', 'Y']

def split_vcf_by_chromosome(directory, vcf_file):
    """Split a VCF file into individual files for each chromosome."""
    dir_name = os.path.basename(directory)
    
    with open(vcf_file, 'r') as file:
        header_lines = []
        chromosome_data = {}
        
        for line in file:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                chromosome = line.split('\t')[0]
                if is_valid_chromosome(chromosome):
                    if chromosome not in chromosome_data:
                        chromosome_data[chromosome] = []
                    chromosome_data[chromosome].append(line)
    
    header = ''.join(header_lines)
    for chromosome, data in chromosome_data.items():
        output_file = os.path.join(directory, f'{dir_name}_chr_{chromosome}.mapped.var.-final.vcf')
        with open(output_file, 'w') as file:
            file.write(header)
            file.write(''.join(data))
        print(f'File for chromosome {chromosome} saved as: {output_file}')

def main():
    current_directory = os.getcwd()
    for directory in glob.glob(os.path.join(current_directory, '[ES]RR*')):
        if os.path.isdir(directory):
            vcf_files = glob.glob(os.path.join(directory, '*var.-final.vcf'))
            if vcf_files:
                vcf_file = vcf_files[0]  # Assuming there's only one such file per directory
                print(f'Processing file: {vcf_file}')
                split_vcf_by_chromosome(directory, vcf_file)
            else:
                print(f'No VCF file found in directory: {directory}')

if __name__ == "__main__":
    main()
