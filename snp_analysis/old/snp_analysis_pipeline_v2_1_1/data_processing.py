import os
import shutil
import re
from snp_analysis_pipeline_v2_1.logging_module import log_message
from snp_analysis_pipeline_v2_1.command_execution import run_command

def detect_read_type(accession_number):
    """Determine read type based on available files."""
    # List all FASTQ files in the accession directory
    fastq_files = [f for f in os.listdir(accession_number) if re.match(f"{accession_number}_\\d+\\.fastq$", f)]
    
    # Check for _trimmed.fq file for trimmed data
    trimmed_file = f"{accession_number}_trimmed.fq"
    if trimmed_file in os.listdir(accession_number):
        fastq_files.append(trimmed_file)
    
    if not fastq_files:
        log_message(f"No FASTQ files found for {accession_number}.", level="error")
        raise FileNotFoundError(f"No FASTQ files found for {accession_number}.")
    
    # Check for paired-end read pattern (e.g., _1.fastq and _2.fastq) or single trimmed file
    paired_files = [f for f in fastq_files if re.match(f"{accession_number}_[12]\\.fastq$", f)]
    
    if len(paired_files) == 2:
        read_type = '2'
        log_message(f"Detected paired-end reads for {accession_number}.", level="info")
    elif trimmed_file in fastq_files:
        read_type = '1'
        log_message(f"Detected single-end trimmed reads for {accession_number}.", level="info")
    else:
        read_type = '1'
        log_message(f"Detected single-end reads for {accession_number}.", level="info")
    
    return read_type

def download_and_convert(accession_number, threads):
    """Download and convert SRA files using fastq-dump."""
    # Skip downloading if trimmed files are already present
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq"):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping download and conversion.", level="info")
        return

    # Ensure the target directory exists
    os.makedirs(accession_number, exist_ok=True)

    log_message(f"\nDownloading and converting {accession_number} with fastq-dump...", level="info")

    run_command(f"fastq-dump --split-files --outdir {accession_number} {accession_number}")

    # Move FASTQ files to the accession directory if they exist outside
    for file in os.listdir('.'):
        if file.startswith(f"{accession_number}_") and file.endswith(".fastq"):
            shutil.move(file, os.path.join(accession_number, file))

    # Verify that files are in place
    fastq_files = [f for f in os.listdir(accession_number) if f.endswith(".fastq")]

    if not fastq_files:
        log_message(f"Failed to find FASTQ files after conversion for {accession_number}", level="error")
        raise FileNotFoundError(f"No FASTQ files found in {accession_number} directory.")


def trim_reads(accession_number, read_type, fastp_path):
    """Trim reads using fastp."""
    # Skip trimming if already trimmed files exist
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq"):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping trimming.", level="info")
        return

    log_message(f"\nTrimming {accession_number} with fastp...", level="info")

    if read_type == '1':
        # For single-end reads, collect all FASTQ files
        fastq_files = [os.path.join(accession_number, f) for f in os.listdir(accession_number) if re.match(f"{accession_number}_\\d+\\.fastq$", f)]
        if not fastq_files:
            log_message(f"No input FASTQ files found for single-end trimming of {accession_number}.", level="error")
            raise FileNotFoundError(f"No input FASTQ files found for {accession_number}.")
        
        # Concatenate all FASTQ files into one
        concatenated_fastq = f"{accession_number}/{accession_number}_combined.fastq"
        with open(concatenated_fastq, 'wb') as wfd:
            for f in fastq_files:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
        
        # Run fastp on the concatenated file
        trim_command = f"{fastp_path} -i {concatenated_fastq} -o {accession_number}/{accession_number}_trimmed.fq --thread=4"
        run_command(trim_command)
    elif read_type == '2':
        input_file_1 = f"{accession_number}/{accession_number}_1.fastq"
        input_file_2 = f"{accession_number}/{accession_number}_2.fastq"
        if not os.path.isfile(input_file_1) or not os.path.isfile(input_file_2):
            log_message(f"One or both input files {input_file_1}, {input_file_2} not found for paired-end trimming.", level="error")
            raise FileNotFoundError(f"Input files for paired-end reads not found.")
        trim_command = f"{fastp_path} -i {input_file_1} -I {input_file_2} -o {accession_number}/{accession_number}_1_trimmed.fq -O {accession_number}/{accession_number}_2_trimmed.fq --thread=4"
        run_command(trim_command)

    shutil.move("fastp.html", f"{accession_number}/fastp.html")
    shutil.move("fastp.json", f"{accession_number}/fastp.json")