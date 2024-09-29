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
    log_message("Would you like to download a normal combined VCF file (yes or no)? ", "blue")
    normal_file = input().strip().lower()
    
    if normal_file == 'yes':
        if chromosome in link_list:
            file_id = link_list[chromosome].split('/d/')[1].split('/view')[0]
            url = f"https://drive.google.com/uc?export=download&id={file_id}"
            response = requests.get(url)
            filename = f"ch{chromosome}_normal_comb.vcf"
            output_path = output_directory / filename
            
            if response.status_code == 200:
                with open(output_path, "wb") as f:
                    f.write(response.content)
                log_message(f"{filename} downloaded successfully!", "success")
            else:
                log_message(f"Failed to download the file. Status code: {response.status_code}", "error")
        else:
            log_message("Chromosome number not found in link list.", "error")
    elif normal_file == 'no':
        log_message("Download skipped.", "info")
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

        # Download the normal VCF file and place it in the parent directory
        link_list = {
            "1": "https://drive.google.com/file/d/1-09an3LHEJYnf-0s4S70ATis81dZPV1G/view?usp=sharing",
            "2": "https://drive.google.com/file/d/1xP8ELyW5Kr7fay0NulbgiPQOAx5Eq9E3/view?usp=sharing",
            "3": "https://drive.google.com/file/d/1-JNXwA4VvwcSOoulcmRmkLq8YueltyQ/view?usp=sharing",
            "4": "https://drive.google.com/file/d/1n2T5xzMOqOobm5w97uthlvW5d7e3nV9M/view?usp=sharing",
            "5": "https://drive.google.com/file/d/1s8OPif3dKAesCgoOTIQGdAr9CdcUhe2q/view?usp=sharing",
            "6": "https://drive.google.com/file/d/1IzA-rUe8ELPyFtCePiCQgQeWDnOEPJxN/view?usp=sharing",
            "7": "https://drive.google.com/file/d/16eGfwwe1yjvyY6TO-hWj1Rx7YMvXuSbO/view?usp=sharing",
            "8": "https://drive.google.com/file/d/1uSqwMSvQpyBxPwiHhgUcNVkL7YdkcOHg/view?usp=sharing",
            "9": "https://drive.google.com/file/d/1chh47ez6Vqe1W0bdoKqxPHO8iOgMMiJe/view?usp=sharing",
            "10": "https://drive.google.com/file/d/12NXqBfLFJZwifp7_LLAzZFsfpE60czgN/view?usp=sharing",
            "11": "https://drive.google.com/file/d/1nsicoTeVg3t4AL592QGLa0QCiOmrhZz4/view?usp=sharing",
            "12": "https://drive.google.com/file/d/1hhX-NaqjzkpsA9cvFMj4F-1cEPCl6JmK/view?usp=sharing",
            "13": "https://drive.google.com/file/d/112aOwZiP4kOJhrBsBockjM0q2i1a1lx-/view?usp=sharing",
            "14": "https://drive.google.com/file/d/1vS_Fn_uVQvCATA402LC22QMnd5gi95pP/view?usp=sharing",
            "15": "https://drive.google.com/file/d/1-umaSF-rAFqqSwVeQI6djyW7QKxBm5A7/view?usp=sharing",
            "16": "https://drive.google.com/file/d/1cOFSi6ujbJV_I144je6PnBUmiAZL1-J0/view?usp=share_link",
            "17": "https://drive.google.com/file/d/1yhKZrhhw2YlBMySWa4L-DYaJ4UVYH-y8/view?usp=share_link",
            "18": "https://drive.google.com/file/d/1yvj7_SQE93sJ2kwceFr5ikYsFstOZ36I/view?usp=sharing",
            "19": "https://drive.google.com/file/d/1hRKEmwtw7uBkx8n_5wBWvdY2WqQjuOZz/view?usp=sharing",
            "20": "https://drive.google.com/file/d/15vb8C0tkFCczQ4uq-VWnsIQBIcxgVl9h/view?usp=sharing",
            "21": "https://drive.google.com/file/d/1kV21OWPOQdG_0mlgAIua1bJNhpf9UleN/view?usp=sharing",
            "22": "https://drive.google.com/file/d/1zqOFpGV8HEVBo_uDIqVi0_YVQOCL6LZZ/view?usp=share_link",
            "X": "https://drive.google.com/file/d/1y9ZlrPOMFz6bUctxe93x7YBchYBkqvIL/view?usp=sharing",
            "Y": "https://drive.google.com/file/d/1sNGm-4emFDuNOKnGjWozTtvLoIMXUpUg/view?usp=sharing",
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
