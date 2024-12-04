import os
import shutil
import subprocess
from pathlib import Path
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
        os.makedirs(output_dir, exist_ok=True)
        print(f'{GREEN}Created output directory: {output_dir}{RESET}')

    # Collect BAM files grouped by chromosome
    chromosome_files = {}
    for dir_name in os.listdir():
        if 'RR' in dir_name and os.path.isdir(dir_name):
            for file_name in os.listdir(dir_name):
                if file_name.endswith('.sorted.bam'):
                    # Extract chromosome number from file name (assumes format like mapped_17.sorted.bam)
                    chromosome = file_name.split('_')[-2]
                    input_bam_path = Path(dir_name) / file_name
                    if chromosome not in chromosome_files:
                        chromosome_files[chromosome] = []
                    chromosome_files[chromosome].append(input_bam_path)

    # Process BAM files chromosome by chromosome
    for chromosome, bam_files in sorted(chromosome_files.items()):  # Sort by chromosome
        print(f'{BLUE}Processing chromosome: {chromosome}{RESET}')
        for input_bam_path in bam_files:
            file_name = input_bam_path.name
            output_csv_file = output_dir / f"{file_name[:-4]}_counts.csv"

            # Skip if output CSV already exists and is not empty
            if output_csv_file.exists() and output_csv_file.stat().st_size > 0:
                print(f'{MAGENTA}Skipping {file_name} as counts CSV file already exists and is not empty.{RESET}')
                continue

            # Check if the BAM file is valid
            print(f'{MAGENTA}Checking BAM file: {input_bam_path}{RESET}')
            if subprocess.call(['samtools', 'quickcheck', str(input_bam_path)]) != 0:
                print(f'{RED}{input_bam_path} is not a valid BAM file{RESET}')
                continue

            # Indexing the BAM file
            if not (input_bam_path.with_suffix('.bai').exists()):
                print(f'{MAGENTA}Indexing BAM file: {input_bam_path}{RESET}')
                subprocess.call(['samtools', 'index', str(input_bam_path)])

            # Generate count matrix for the mapped BAM file using HTSeq
            print(f'{MAGENTA}Generating count matrix for {input_bam_path}...{RESET}')
            cmd = f'{htseq_count_path} -f bam -s no -i gene_id {input_bam_path} {annotation_file}'
            try:
                with open(output_csv_file, 'w', newline='') as f_out:
                    subprocess.run(cmd.split(), stdout=f_out, check=True)
            except subprocess.CalledProcessError as e:
                print(f'{RED}Error while running htseq-count: {e}{RESET}')

if __name__ == "__main__":
    annotation_file_path = input(f'{MAGENTA}Enter the path to the annotation file in GTF format:{RESET} ')
    output_dir_path = input(f'{MAGENTA}\nEnter the path to the output directory:{RESET} ')
    main(annotation_file_path, output_dir_path)
