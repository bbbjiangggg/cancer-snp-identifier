#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import shutil
import requests
import pysam
import re  # Added for regex matching

# Define color-coded logging
def log_message(message, level="normal", normal_part=""):
    color_map = {
        "info": "\033[95m",       # Magenta
        "warning": "\033[93m",    # Yellow
        "error": "\033[91m",      # Red
        "success": "\033[92m",    # Green
        "normal": "\033[0m",      # Reset to default
        "blue": "\033[94m"        # Blue
    }
    color = color_map.get(level, "\033[0m")
    reset = "\033[0m"
    print(f"{color}{message}{reset}{normal_part}")

def get_parent_directory():
    return Path.cwd()

def list_subdirectories(parent_directory):
    return [d for d in parent_directory.iterdir() if d.is_dir()]

def check_subdirectory_for_vcf(subdir):
    return any(file.suffix in [".vcf", ".vcf.gz"] and "var.-final" in file.stem for file in subdir.iterdir())

def prompt_remove_directory(subdir):
    remove = input(log_message("Remove? (yes/no): ", "blue"))
    if remove.lower() in ("yes", "y"):
        shutil.rmtree(subdir)
        log_message(f"Removed subdirectory {subdir}.", "success")
    else:
        log_message("Subdirectory not removed.", "info")

def process_subdirectories(subdirectories):
    for subdir in subdirectories:
        if subdir.name.startswith(("SRR", "ERR")):
            if check_subdirectory_for_vcf(subdir):
                log_message(f"Verified: {subdir.name}", "info")
            else:
                log_message(f"No file ending with 'var.-final.vcf' found in subdirectory: {subdir.name}", "warning")
                prompt_remove_directory(subdir)

def create_directory(directory):
    if directory.exists():
        log_message("Directory already exists, removing directory.", "error")
        shutil.rmtree(directory)
    directory.mkdir()
    log_message(f"New directory has been created: {directory.name}", "success")

def copy_vcf_files(src_dir, dest_dir, chromosome):
    pattern = re.compile(rf'.*_{chromosome}\.var\.-final\.vcf(\.gz)?$')
    for subdir in src_dir.iterdir():
        if 'RR' in subdir.name:
            for file in subdir.iterdir():
                if file.suffix in [".vcf", ".vcf.gz"] and pattern.match(file.name):
                    shutil.copy2(str(file), dest_dir)
                    log_message(f"Copied {file.name} to {dest_dir}", "info")

def compress_vcf_files(directory):
    for file in directory.iterdir():
        if file.suffix == ".vcf":
            compressed_file = file.with_suffix(".vcf.gz")
            log_message(f"Compressing {file.name} with bgzip...", "info")
            try:
                subprocess.run(["bgzip", "-f", str(file)], check=True)
                log_message(f"Compressed {file.name} to {compressed_file.name}", "success")
            except subprocess.CalledProcessError as e:
                log_message(f"Error compressing {file.name}: {str(e)}", "error")

def index_individual_vcf_files(directory):
    for file in directory.iterdir():
        if file.suffix == ".gz":  # Only index the compressed VCF files
            log_message(f"Indexing {file.name} with tabix...", "info")
            try:
                subprocess.run(["tabix", "-p", "vcf", str(file)], check=True)
                log_message(f"Indexed {file.name}", "success")
            except subprocess.CalledProcessError as e:
                log_message(f"Error indexing {file.name}: {str(e)}", "error")

def create_backup(directory):
    copy_name = f"copy_{directory.name}"
    copy_dir = directory.parent / copy_name
    if copy_dir.exists():
        log_message("Backup directory already exists, removing directory.", "error")
        shutil.rmtree(copy_dir)
    shutil.copytree(directory, copy_dir)
    log_message(f"A copy of {directory.name} has been created.", "success")

def merge_vcf_files(vcf_files, output_file):
    merge_command = ["bcftools", "merge", "-o", output_file, "-O", "z"] + vcf_files
    try:
        subprocess.run(merge_command, check=True)
        log_message(f"Combined VCF file created: {output_file}", "success")
    except subprocess.CalledProcessError as e:
        log_message(f"Error during merging: {str(e)}", "error")

def index_vcf_file(vcf_file):
    try:
        subprocess.run(["tabix", "-p", "vcf", vcf_file], check=True)
        log_message(f"VCF file indexed: {vcf_file}", "success")
    except subprocess.CalledProcessError as e:
        log_message(f"Error during indexing: {str(e)}", "error")

def download_normal_vcf(link_list, output_directory, chromosome, cancer):
    log_message("Would you like to copy a normal combined VCF file (yes or no)? ", "blue")
    normal_file = input().strip().lower()
    
    if normal_file == 'yes':
        if chromosome in link_list:
            filename = f"ch{chromosome}_normal_comb.vcf"
            source_path = link_list[chromosome] / filename
            output_path = output_directory / filename
            
            if source_path.exists():
                shutil.copy2(str(source_path), str(output_path))
                log_message(f"{filename} copied successfully!", "success")
            else:
                log_message(f"File {filename} not found in the source directory.", "error")
        else:
            log_message("Chromosome number not found in link list.", "error")
    elif normal_file == 'no':
        log_message("Copy skipped.", "info")
    else:
        log_message("Invalid input. Please enter 'yes' or 'no'.", "error")

def generate_binary_string_vcf(input_vcf, output_file):
    vcf_in = pysam.VariantFile(input_vcf)  # Open the VCF file with pysam
    with open(output_file, 'w') as out_f:
        for record in vcf_in:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0] if record.alts else "."
            
            binary_string = ""
            for sample in record.samples.values():
                gt = sample['GT']
                if gt == (1, 1) or gt == (0, 1):
                    binary_string += '1'
                else:
                    binary_string += '0'
                    
            out_f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{binary_string}\n")
    log_message(f"Binary string VCF-like file saved to {output_file}", "success")

def move_files_to_directory(parent_directory, chromosome, cancer, files):
    output_dir = parent_directory / f"ch{chromosome}_{cancer}_output"
    create_directory(output_dir)

    possible_source_dirs = [
        parent_directory,
        parent_directory / f"ch{chromosome}_{cancer}_vcf"
    ]

    for file in files:
        found = False
        for source_dir in possible_source_dirs:
            file_path = source_dir / file
            if file_path.exists():
                shutil.move(str(file_path), str(output_dir / file_path.name))
                log_message(f"Moved {file_path.name} from {source_dir} to {output_dir}", "success")
                found = True
                break
        if not found:
            log_message(f"File {file} not found in any source directory, skipping move.", "warning")

def main():
    parent_directory = get_parent_directory()
    log_message(f"Parent directory path: {parent_directory} ", "info")

    subdirectories = list_subdirectories(parent_directory)
    log_message("Checking subdirectories for var.-final.vcf...", "info")
    process_subdirectories(subdirectories)

    log_message("Enter the chromosome number of this analysis (1-22, X, or Y): ", "blue")
    chromosome = input().strip()

    log_message("Enter an abbreviation for the cancer type (e.g. pnca): ", "blue")
    cancer = input().strip()

    directory_name = f'ch{chromosome}_{cancer}_vcf'
    isecdir = parent_directory / directory_name

    create_directory(isecdir)

    log_message("Press enter to continue...", "blue")
    input()

    log_message("Processing files, please wait...", "info")

    copy_vcf_files(parent_directory, isecdir, chromosome)
    compress_vcf_files(isecdir)
    index_individual_vcf_files(isecdir)
    create_backup(isecdir)

    os.chdir(isecdir)
    vcf_files = sorted([str(file) for file in isecdir.glob("*.vcf.gz")])

    if vcf_files:
        output_file = f"ch{chromosome}_{cancer}_comb.vcf.gz"
        merge_vcf_files(vcf_files, output_file)
        index_vcf_file(output_file)

        # Generate the combined VCF (previously binary string VCF-like file)
        combined_vcf_file = parent_directory / f"ch{chromosome}_{cancer}_comb.vcf"
        generate_binary_string_vcf(output_file, combined_vcf_file)

        # Define the source directory for normal VCF files
        source_directory = Path("/usr/local/bin/normal_vcf")
        link_list = {
            "1": source_directory,
            "2": source_directory,
            "3": source_directory,
            "4": source_directory,
            "5": source_directory,
            "6": source_directory,
            "7": source_directory,
            "8": source_directory,
            "9": source_directory,
            "10": source_directory,
            "11": source_directory,
            "12": source_directory,
            "13": source_directory,
            "14": source_directory,
            "15": source_directory,
            "16": source_directory,
            "17": source_directory,
            "18": source_directory,
            "19": source_directory,
            "20": source_directory,
            "21": source_directory,
            "22": source_directory,
            "X": source_directory,
            "Y": source_directory,
        }
        download_normal_vcf(link_list, parent_directory, chromosome, cancer)

        # Move the generated files into the output directory
        files_to_move = [
            f"ch{chromosome}_normal_comb.vcf",
            f"ch{chromosome}_{cancer}_comb.vcf",
            f"ch{chromosome}_{cancer}_comb.vcf.gz",
            f"ch{chromosome}_{cancer}_comb.vcf.gz.tbi"
        ]

        move_files_to_directory(parent_directory, chromosome, cancer, files_to_move)
    else:
        log_message("No VCF files found for merging.", "warning")

if __name__ == "__main__":
    main()