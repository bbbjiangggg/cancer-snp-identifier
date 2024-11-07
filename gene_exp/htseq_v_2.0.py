
import os
import shutil
import subprocess
from pathlib import Path
import csv
import gzip

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def detect_htseq_path():
    htseq_path = shutil.which('htseq-count')
    if htseq_path is None:
        print(f'{RED}htseq-count not found in PATH.{RESET}')
        print(f'{MAGENTA}Please run the command "which htseq-count" to find the path to htseq-count and enter it below.{RESET}')
        htseq_path = input(f'{MAGENTA}Enter the full path to htseq-count:{RESET} ')
    return htseq_path

def main(annotation_file_path: str, output_dir_path: str):
    annotation_file = Path(annotation_file_path)
    output_dir = Path(output_dir_path)
    htseq_count_path = detect_htseq_path()

    # Decompress the GTF file if it is in .gz format
    if annotation_file_path.endswith('.gz'):
        decompressed_file_path = annotation_file_path[:-3]  # Remove .gz extension
        print(f'{MAGENTA}Decompressing {annotation_file_path}...{RESET}')
        with gzip.open(annotation_file_path, 'rt') as gz_file, open(decompressed_file_path, 'w') as output_file:
            output_file.write(gz_file.read())
        annotation_file_path = decompressed_file_path

    if not annotation_file.is_file():
        print(f'{RED}{annotation_file} is not a valid file path{RESET}')
        return

    if not output_dir.is_dir():
        print(f'{RED}{output_dir} is not a valid directory path{RESET}')
        return

    # Loop through all directories that have 'RR' in the name
    for dir_name in os.listdir():
        if 'RR' in dir_name and os.path.isdir(dir_name):
            print(f'{BLUE}Processing directory: {dir_name}{RESET}')
            # Search for the mapped and sorted BAM file in the directory
            for file_name in os.listdir(dir_name):
                if file_name.endswith('.sorted.bam'):
                    print(f'{GREEN}Found BAM file: {file_name}{RESET}')
                    output_csv_file = output_dir / f"{file_name[:-4]}_counts.csv"
                    
                    # Skip if output CSV already exists and is not empty
                    if output_csv_file.exists() and output_csv_file.stat().st_size > 0:
                        print(f'{MAGENTA}Skipping {file_name} as counts CSV file already exists and is not empty.{RESET}')
                        continue

                    output_bam_path = output_dir / file_name
                    input_bam_path = Path(dir_name) / file_name

                    # Check if the BAM file is valid
                    print(f'{MAGENTA}Checking BAM file: {input_bam_path}{RESET}')
                    if subprocess.call(['samtools', 'quickcheck', str(input_bam_path)]) != 0:
                        print(f'{RED}{input_bam_path} is not a valid BAM file{RESET}')
                        continue

                    # Check if the BAM file is sorted (already ensured by the naming convention)
                    header = subprocess.check_output(['samtools', 'view', '-H', str(input_bam_path)])
                    if b'SO:coordinate' not in header:
                        print(f'{RED}{input_bam_path} is not sorted by genomic coordinates{RESET}')
                        # Sort the BAM file
                        sort_bam_path = output_dir / f"{file_name[:-4]}_sorted.bam"
                        print(f'{MAGENTA}Sorting {input_bam_path} to {sort_bam_path}...{RESET}')
                        sort_cmd = f'samtools sort -o {sort_bam_path} {input_bam_path}'
                        subprocess.call(sort_cmd, shell=True)
                        input_bam_path = sort_bam_path

                    # Indexing the BAM file
                    print(f'{MAGENTA}Indexing BAM file: {output_bam_path}{RESET}')
                    subprocess.call(['samtools', 'index', str(output_bam_path)])

                    # Generate count matrix for the mapped BAM file using HTSeq
                    print(f'{MAGENTA}Generating count matrix for {output_bam_path}...{RESET}')
                    cmd = f'{htseq_count_path} -f bam -s no -i gene_id {output_bam_path} {annotation_file}'
                    
                    # Run htseq-count and directly save output as CSV
                    with open(output_csv_file, 'w', newline='') as f_out:
                        subprocess.run(cmd.split(), stdout=f_out)

if __name__ == "__main__":
    annotation_file_path = input(f'{MAGENTA}Enter the path to the annotation file in GTF format:{RESET} ')
    output_dir_path = input(f'{MAGENTA}\nEnter the path to the output directory:{RESET} ')
    main(annotation_file_path, output_dir_path)
