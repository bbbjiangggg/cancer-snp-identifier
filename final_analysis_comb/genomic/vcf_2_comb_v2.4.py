#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
import shutil
import requests
import pysam
import re

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
        log_message(f"Directory {directory} already exists, removing directory.", "error")
        shutil.rmtree(directory)
    directory.mkdir()
    log_message(f"New directory has been created: {directory.name}", "success")

def copy_vcf_files(src_dir, dest_dir, chromosome):
    if chromosome.lower() == 'genomic':
        pattern = re.compile(r'^SRR\d+_mapped\.var\.-final\.vcf(\.gz)?$')
    else:
        pattern = re.compile(rf'.*_{chromosome}\.var\.-final\.vcf(\.gz)?$')
    copied_files = []
    for subdir in src_dir.iterdir():
        if 'RR' in subdir.name:
            for file in subdir.iterdir():
                if file.suffix in [".vcf", ".vcf.gz"] and pattern.match(file.name):
                    file_size = file.stat().st_size
                    if file_size == 0:
                        log_message(f"File {file.name} is zero bytes and will be deleted.", "warning")
                        file.unlink()  # Delete the zero-byte file
                        continue
                    shutil.copy2(str(file), dest_dir)
                    log_message(f"Copied {file.name} to {dest_dir}", "info")
                    copied_files.append(file.name)
    if not copied_files:
        log_message("No VCF files were copied. Please check your input files and chromosome selection.", "error")
        exit(1)

def compress_vcf_files(directory):
    for file in directory.iterdir():
        if file.suffix == ".vcf":
            file_size = file.stat().st_size
            if file_size == 0:
                log_message(f"File {file.name} is zero bytes and will be deleted.", "warning")
                file.unlink()  # Delete the zero-byte file
                continue
            compressed_file = file.with_suffix(".vcf.gz")
            log_message(f"Compressing {file.name} with bgzip...", "info")
            try:
                subprocess.run(["bgzip", "-f", str(file)], check=True)
                log_message(f"Compressed {file.name} to {compressed_file.name}", "success")
            except subprocess.CalledProcessError as e:
                log_message(f"Error compressing {file.name}: {str(e)}", "error")
                exit(1)
        elif file.suffix == ".gz":
            # Verify that the file is properly compressed
            result = subprocess.run(["bgzip", "-t", str(file)])
            if result.returncode != 0:
                log_message(f"File {file.name} is not properly compressed. Re-compressing...", "warning")
                try:
                    subprocess.run(["bgzip", "-f", str(file.with_suffix(""))], check=True)
                    log_message(f"Re-compressed {file.name}", "success")
                except subprocess.CalledProcessError as e:
                    log_message(f"Error re-compressing {file.name}: {str(e)}", "error")
                    exit(1)

def index_individual_vcf_files(directory):
    for file in directory.iterdir():
        if file.suffix == ".gz":  # Only index the compressed VCF files
            log_message(f"Indexing {file.name} with tabix...", "info")
            try:
                subprocess.run(["tabix", "-p", "vcf", str(file)], check=True)
                log_message(f"Indexed {file.name}", "success")
            except subprocess.CalledProcessError as e:
                log_message(f"Error indexing {file.name}: {str(e)}", "error")
                exit(1)

def create_backup(directory):
    copy_name = f"copy_{directory.name}"
    copy_dir = directory.parent / copy_name
    if copy_dir.exists():
        log_message("Backup directory already exists, removing directory.", "error")
        shutil.rmtree(copy_dir)
    shutil.copytree(directory, copy_dir)
    log_message(f"A copy of {directory.name} has been created.", "success")

def assign_unique_sample_names(directory):
    for idx, file in enumerate(directory.glob("*.vcf.gz")):
        sample_name = f"Sample_{idx+1}"
        log_message(f"Assigning unique sample name '{sample_name}' to {file.name}", "info")
        try:
            # Extract existing sample names
            original_sample_names = subprocess.check_output(
                ["bcftools", "query", "-l", str(file)],
                universal_newlines=True
            ).strip()
            # Create a mapping from original to new sample name
            with open("samples.txt", "w") as f:
                f.write(f"{original_sample_names}\t{sample_name}\n")
            # Reheader the VCF file
            subprocess.run(
                f"bcftools reheader -s samples.txt -o {file}.tmp {file}",
                shell=True,
                check=True
            )
            os.replace(f"{file}.tmp", file)
            os.remove("samples.txt")
        except subprocess.CalledProcessError as e:
            log_message(f"Error assigning unique sample name to {file.name}: {str(e)}", "error")
            # Move the problematic file to a separate directory for review
            problematic_dir = directory / "problematic_files"
            problematic_dir.mkdir(exist_ok=True)
            shutil.move(str(file), problematic_dir)
            log_message(f"Moved problematic file {file.name} to {problematic_dir}", "warning")

def merge_vcf_files(vcf_files, output_file):
    log_message(f"Merging {len(vcf_files)} VCF files into {output_file}...", "info")
    merge_command = [
        "bcftools", "merge",
        "--force-samples",
        "--merge", "all",
        "-O", "z",
        "-o", output_file
    ] + vcf_files
    try:
        subprocess.run(merge_command, check=True)
        log_message(f"Combined VCF file created: {output_file}", "success")
    except subprocess.CalledProcessError as e:
        log_message(f"Error during merging: {str(e)}", "error")
        exit(1)

def index_vcf_file(vcf_file):
    log_message(f"Indexing combined VCF file {vcf_file}...", "info")
    try:
        subprocess.run(["tabix", "-p", "vcf", vcf_file], check=True)
        log_message(f"VCF file indexed: {vcf_file}", "success")
    except subprocess.CalledProcessError as e:
        log_message(f"Error during indexing: {str(e)}", "error")
        exit(1)

def download_normal_vcf(link_list, output_directory, chromosome, cancer):
    # Skip downloading the normal VCF file when 'genomic' is selected
    if chromosome.lower() == 'genomic':
        return
    log_message("Would you like to download a normal combined VCF file (yes or no)? ", "blue")
    normal_file = input().strip().lower()
    
    if normal_file == 'yes':
        if chromosome in link_list:
            file_id = link_list[chromosome].split('/d/')[1].split('/view')[0]
            url = f"https://drive.google.com/uc?export=download&id={file_id}"
            log_message(f"Downloading normal VCF file for chromosome {chromosome}...", "info")
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
    if not os.path.exists(input_vcf):
        log_message(f"Combined VCF file {input_vcf} does not exist.", "error")
        exit(1)
    vcf_in = pysam.VariantFile(input_vcf)  # Open the VCF file with pysam
    with open(output_file, 'w') as out_f:
        for record in vcf_in:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = ','.join(record.alts) if record.alts else "."
            
            binary_string = ""
            for sample in vcf_in.header.samples:
                sample_data = record.samples[sample]
                gt = sample_data.get('GT')
                if gt in [(0, 1), (1, 1), (1,)]:
                    binary_string += '1'
                elif gt == (None, None) or gt == ('.', '.'):
                    binary_string += '0'  # Treat missing data as '0'
                else:
                    binary_string += '0'
                    
            out_f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{binary_string}\n")
    log_message(f"Binary string VCF-like file saved to {output_file}", "success")

def move_files_to_directory(parent_directory, chromosome_str, cancer, files):
    output_dir = parent_directory / f"{chromosome_str}_{cancer}_output"
    create_directory(output_dir)

    possible_source_dirs = [
        parent_directory,
        parent_directory / f"{chromosome_str}_{cancer}_vcf"
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

    log_message("Enter the chromosome number of this analysis (1-22, X, Y, or 'genomic'): ", "blue")
    chromosome = input().strip()

    log_message("Enter an abbreviation for the cancer type (e.g. pnca): ", "blue")
    cancer = input().strip()

    if chromosome.lower() == 'genomic':
        chromosome_str = 'genomic'
    else:
        chromosome_str = f"ch{chromosome}"

    directory_name = f'{chromosome_str}_{cancer}_vcf'
    isecdir = parent_directory / directory_name

    create_directory(isecdir)

    log_message("Press enter to continue...", "blue")
    input()

    log_message("Processing files, please wait...", "info")

    copy_vcf_files(parent_directory, isecdir, chromosome)
    compress_vcf_files(isecdir)
    index_individual_vcf_files(isecdir)
    create_backup(isecdir)

    # Assign unique sample names before merging
    assign_unique_sample_names(isecdir)

    os.chdir(isecdir)
    vcf_files = sorted([str(file) for file in isecdir.glob("*.vcf.gz")])

    if vcf_files:
        output_file = f"{chromosome_str}_{cancer}_comb.vcf.gz"
        merge_vcf_files(vcf_files, output_file)

        if os.path.exists(output_file):
            index_vcf_file(output_file)
        else:
            log_message(f"Merged VCF file {output_file} was not created. Merging failed.", "error")
            exit(1)

        # Generate the combined VCF (previously binary string VCF-like file)
        combined_vcf_file = parent_directory / f"{chromosome_str}_{cancer}_comb.vcf"
        generate_binary_string_vcf(output_file, combined_vcf_file)

        # Skip downloading the normal VCF file when 'genomic' is selected
        if chromosome.lower() != 'genomic':
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
                f"{chromosome_str}_normal_comb.vcf",
                f"{chromosome_str}_{cancer}_comb.vcf",
                f"{chromosome_str}_{cancer}_comb.vcf.gz",
                f"{chromosome_str}_{cancer}_comb.vcf.gz.tbi"
            ]
        else:
            # For 'genomic', move the combined VCF files
            files_to_move = [
                f"{chromosome_str}_{cancer}_comb.vcf",
                f"{chromosome_str}_{cancer}_comb.vcf.gz",
                f"{chromosome_str}_{cancer}_comb.vcf.gz.tbi"
            ]

        move_files_to_directory(parent_directory, chromosome_str, cancer, files_to_move)
    else:
        log_message("No VCF files found for merging.", "warning")
        exit(1)

if __name__ == "__main__":
    main()
