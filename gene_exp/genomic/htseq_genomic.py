#!/usr/bin/env python3

import os
import shutil
import subprocess
from pathlib import Path
import csv



# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def main(annotation_file_path: str, output_dir_path: str):
    annotation_file = Path(annotation_file_path)
    output_dir = Path(output_dir_path)

    if not annotation_file.is_file():
        print(f'{RED}{annotation_file} is not a valid file path{RESET}')
        return

    if not output_dir.is_dir():
        print(f'{RED}{output_dir} is not a valid directory path{RESET}')
        return

    # Loop through all directories that have 'RR' in the name
    for dir_name in os.listdir():
        if 'RR' in dir_name and os.path.isdir(dir_name):
            # Search for the mapped BAM file in the directory
            for file_name in os.listdir(dir_name):
                if file_name.endswith('_mapped_hg38.sorted.bam'):
                    # Move the mapped BAM file to the output directory
                    output_bam_path = output_dir / file_name
                    input_bam_path = Path(dir_name) / file_name

                    # Check if the BAM file is valid
                    if subprocess.call(['samtools', 'quickcheck', str(input_bam_path)]) != 0:
                        print(f'{RED}{input_bam_path} is not a valid BAM file{RESET}')
                        output_count_file = None
                        continue

                    # Check if the BAM file is sorted
                    header = subprocess.check_output(['samtools', 'view', '-H', str(input_bam_path)])
                    if b'SO:coordinate' not in header:
                        print(f'{RED}{input_bam_path} is not sorted by genomic coordinates{RESET}')
                        # Sort the BAM file
                        sort_bam_path = output_dir / f"{file_name[:-4]}_mapped_hg38.sorted.bam"
                        print(f'{MAGENTA}Sorting {input_bam_path} to {sort_bam_path}...{RESET}')
                        sort_cmd = f'samtools sort -o {sort_bam_path} {input_bam_path}'
                        subprocess.call(sort_cmd, shell=True)
                        input_bam_path = sort_bam_path

                        # Check for errors in the alignment after sorting
                        flagstat = subprocess.check_output(['samtools', 'flagstat', str(input_bam_path)])
                        if b'error' in flagstat:
                            print(f'{RED}Error in alignment of {input_bam_path}{RESET}')
                            output_count_file = None
                            continue

                    else:
                        # Check for errors in the alignment before copying the file
                        flagstat = subprocess.check_output(['samtools', 'flagstat', str(input_bam_path)])
                        if b'error' in flagstat:
                            print(f'{RED}Error in alignment of {input_bam_path}{RESET}')
                            output_count_file = None
                            continue

                    output_count_file = output_dir / f"{file_name[:-4]}_counts.txt"
                    shutil.copy(str(input_bam_path), str(output_bam_path))

                    

                    if output_count_file:
                            # Generate count matrix for the mapped BAM file using HTSeq
                            # Create index file
                            subprocess.call(['samtools', 'index', str(output_bam_path)])
                            # Generate count matrix as a text file
                            print(f'{MAGENTA}Generating count matrix for {output_bam_path}...{RESET}')
                            cmd = f'/home/crisprmax/.local/bin/htseq-count -f bam -s no -i gene_id {output_bam_path} {annotation_file} > {output_count_file}'



                            # Run the htseq-count command
                            subprocess.call(cmd, shell=True)

                            # Convert the count matrix to a CSV file
                            csv_count_file = output_count_file.with_suffix('.csv')
                            with open(output_count_file, 'r') as f_in, open(csv_count_file, 'w', newline='') as f_out:
                                csv_reader = csv.reader(f_in, delimiter='\t')
                                csv_writer = csv.writer(f_out)
                                for row in csv_reader:
                                    csv_writer.writerow(row)
                                    
                            # Remove the original count file
                            os.remove(str(output_count_file))

                
if __name__ == "__main__":
    annotation_file_path = input(f'{MAGENTA}Enter the path to the annotation file in GTF format:{RESET} ')
    output_dir_path = input(f'{MAGENTA}\nEnter the path to the output directory:{RESET} ')
    main(annotation_file_path, output_dir_path)
