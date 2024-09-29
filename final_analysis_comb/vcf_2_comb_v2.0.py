#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import shutil
import requests
import io

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
    return any(file.suffix == ".vcf" and "var.-final" in file.stem for file in subdir.iterdir())

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

def copy_vcf_files(src_dir, dest_dir):
    for subdir in src_dir.iterdir():
        if 'RR' in subdir.name:
            for file in subdir.iterdir():
                if file.suffix == ".vcf":
                    shutil.copy2(str(file), dest_dir)

def create_backup(directory):
    copy_name = f"copy_{directory.name}"
    copy_dir = directory.parent / copy_name
    if copy_dir.exists():
        log_message("Directory already exists, removing directory.", "error")
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

def download_normal_vcf(link_list):
    while True:
        normal_file = input(log_message("Would you like to download a normal combined VCF file (yes or no)? ", "blue"))
        if normal_file.lower() == 'yes':
            chromosome = input(log_message("Enter the chromosome number of interest (1-22, X or Y): ", "blue"))
            if chromosome in link_list:
                url = link_list[chromosome].replace('/view', '/uc?export=download')
                response = requests.get(url)
                filename = f"ch{chromosome}_norm_comb.vcf"
                with open(filename, "wb") as f:
                    f.write(response.content)
                shutil.copy2(filename, get_parent_directory())
                log_message(f"{filename} downloaded successfully!", "success")
                break
            else:
                log_message("Chromosome number not found in link list.", "error")
        elif normal_file.lower() == 'no':
            break
        else:
            log_message("Invalid input. Please enter 'yes' or 'no'.", "error")

def main():
    parent_directory = get_parent_directory()
    log_message(f"Parent directory path: {parent_directory} ", "info")

    subdirectories = list_subdirectories(parent_directory)
    log_message("Checking subdirectories for var.-final.vcf...", "info")
    process_subdirectories(subdirectories)

    chr = input(log_message("Enter the chromosome number of this analysis (1-22, X, or Y): ", "blue"))
    cancer = input(log_message("Enter an abbreviation for the cancer type (e.g. pnca): ", "blue"))
    directory_name = f'ch{chr}_{cancer}_vcf'
    isecdir = parent_directory / directory_name

    create_directory(isecdir)
    input(log_message("Press enter to continue...", "blue"))

    copy_vcf_files(parent_directory, isecdir)
    create_backup(isecdir)

    os.chdir(isecdir)
    vcf_files = sorted([str(file) for file in isecdir.glob("*.vcf.gz")])

    if vcf_files:
        output_file = f"ch{chr}_{cancer}_comb.vcf.gz"
        merge_vcf_files(vcf_files, output_file)
        index_vcf_file(output_file)
    else:
        log_message("No VCF files found for merging.", "warning")

    link_list = {
        "1": "https://drive.google.com/file/d/1-09an3LHEJYnf-0s4S70ATis81dZPV1G/view?usp=sharing",
        "2": "https://drive.google.com/file/d/1xP8ELyW5Kr7fay0NulbgiPQOAx5Eq9E3/view?usp=sharing",
        "3": "https://drive.google.com/file/d/1-JNXwA4VvwcSOoulcmRmkLq8qYueltyQ/view?usp=sharing",
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
    download_normal_vcf(link_list)

if __name__ == "__main__":
    main()
