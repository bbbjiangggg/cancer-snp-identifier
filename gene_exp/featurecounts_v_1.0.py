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

def detect_featurecounts_path():
    featurecounts_path = shutil.which('featureCounts')
    if featurecounts_path is None:
        print(f'{RED}featureCounts not found in PATH.{RESET}')
        print(f'{MAGENTA}Please run the command "which featureCounts" to find the path to featureCounts and enter it below.{RESET}')
        featurecounts_path = input(f'{MAGENTA}Enter the full path to featureCounts:{RESET} ')
    return featurecounts_path

def main(annotation_file_path: str, output_dir_path: str):
    annotation_file = Path(annotation_file_path)
    output_dir = Path(output_dir_path)
    featurecounts_path = detect_featurecounts_path()

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
                    # Move the mapped BAM file to the output directory
                    output_bam_path = output_dir / file_name
                    input_bam_path = Path(dir_name) / file_name

                    # Check if the BAM file is valid
                    print(f'{MAGENTA}Checking BAM file: {input_bam_path}{RESET}')
                    if subprocess.call(['samtools', 'quickcheck', str(input_bam_path)]) != 0:
                        print(f'{RED}{input_bam_path} is not a valid BAM file{RESET}')
                        output_count_file = None
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

                        # Check for errors in the alignment after sorting
                        print(f'{MAGENTA}Checking alignment for errors: {input_bam_path}{RESET}')
                        flagstat = subprocess.check_output(['samtools', 'flagstat', str(input_bam_path)])
                        if b'error' in flagstat:
                            print(f'{RED}Error in alignment of {input_bam_path}{RESET}')
                            output_count_file = None
                            continue

                    else:
                        # Check for errors in the alignment before copying the file
                        print(f'{MAGENTA}Checking alignment for errors: {input_bam_path}{RESET}')
                        flagstat = subprocess.check_output(['samtools', 'flagstat', str(input_bam_path)])
                        if b'error' in flagstat:
                            print(f'{RED}Error in alignment of {input_bam_path}{RESET}')
                            output_count_file = None
                            continue

                    output_count_file = output_dir / f"{file_name[:-4]}_counts.txt"
                    print(f'{MAGENTA}Copying BAM file to output directory: {output_bam_path}{RESET}')
                    shutil.copy(str(input_bam_path), str(output_bam_path))

                    if output_count_file:
                        # Generate count matrix for the mapped BAM file using featureCounts
                        # Create index file
                        print(f'{MAGENTA}Indexing BAM file: {output_bam_path}{RESET}')
                        subprocess.call(['samtools', 'index', str(output_bam_path)])
                        # Generate count matrix as a text file
                        print(f'{MAGENTA}Generating count matrix for {output_bam_path}...{RESET}')
                        cmd = f'{featurecounts_path} -T 4 -t exon -g gene_id -a {annotation_file} -o {output_count_file} {output_bam_path}'

                        # Run the featureCounts command
                        subprocess.call(cmd, shell=True)

                        # Convert the count matrix to a CSV file
                        csv_count_file = output_count_file.with_suffix('.csv')
                        print(f'{MAGENTA}Converting count matrix to CSV: {csv_count_file}{RESET}')
                        with open(output_count_file, 'r') as f_in, open(csv_count_file, 'w', newline='') as f_out:
                            csv_reader = csv.reader(f_in, delimiter='\t')
                            csv_writer = csv.writer(f_out)
                            for row in csv_reader:
                                csv_writer.writerow(row)

                        # Remove the original count file
                        print(f'{MAGENTA}Removing original count file: {output_count_file}{RESET}')
                        os.remove(str(output_count_file))

if __name__ == "__main__":
    annotation_file_path = input(f'{MAGENTA}Enter the path to the annotation file in GTF format:{RESET} ')
    output_dir_path = input(f'{MAGENTA}\nEnter the path to the output directory:{RESET} ')
    main(annotation_file_path, output_dir_path)
