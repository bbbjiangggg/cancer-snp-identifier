import os
import shutil
from snp_analysis_pipeline_v2_1.logging_module import log_message
from snp_analysis_pipeline_v2_1.command_execution import run_command

def detect_read_type(accession_number):
    """Automatically detect if reads are single-end or paired-end."""
    # Check if trimmed files exist and skip FASTQ detection
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq")):
        log_message(f"Trimmed files detected for {accession_number}. Skipping FASTQ detection.", level="info")
        if os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq"):
            return '2'
        else:
            return '1'

    # If no trimmed files are found, check for FASTQ files
    if os.path.isfile(f"{accession_number}/{accession_number}_1.fastq") and os.path.isfile(f"{accession_number}/{accession_number}_2.fastq"):
        log_message(f"Detected paired-end reads for {accession_number}.", level="info")
        return '2'
    elif os.path.isfile(f"{accession_number}/{accession_number}.fastq"):
        log_message(f"Detected single-end reads for {accession_number}.", level="info")
        return '1'
    else:
        log_message(f"No valid FASTQ files found for {accession_number}.", level="error")
        raise FileNotFoundError(f"No FASTQ files found for {accession_number}.")

def prefetch_and_convert(accession_number, threads):
    """Run prefetch and fasterq-dump to download and convert SRA files."""
    # Ask user if they want to use fastq-dump instead of prefetch and fasterq-dump
    use_fastq_dump = input("Would you like to use fastq-dump instead of prefetch and fasterq-dump? (yes/y or no/n): ").strip().lower()
    
    # Skip downloading if trimmed files are already present
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq")):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping download and conversion.", level="info")
        return

    if use_fastq_dump in ['yes', 'y']:
        log_message(f"\nDownloading and converting {accession_number} with fastq-dump...", level="info")
        if detect_read_type(accession_number) == '2':
            run_command(f"fastq-dump --split-files {accession_number} --outdir {accession_number}")
        else:
            run_command(f"fastq-dump {accession_number} --outdir {accession_number}")
    else:
        log_message(f"\nDownloading sequence number {accession_number} with prefetch...", level="info")
        run_command(f"prefetch {accession_number}")

        log_message(f"\nConverting {accession_number} with fasterq-dump...", level="info")
        run_command(f"fasterq-dump {accession_number} --threads {threads} --progress")

    # Move fastq files to accession directory if not already there
    if os.path.isfile(f"{accession_number}.fastq"):
        shutil.move(f"{accession_number}.fastq", f"{accession_number}/{accession_number}.fastq")
    if os.path.isfile(f"{accession_number}_1.fastq"):
        shutil.move(f"{accession_number}_1.fastq", f"{accession_number}/{accession_number}_1.fastq")
    if os.path.isfile(f"{accession_number}_2.fastq"):
        shutil.move(f"{accession_number}_2.fastq", f"{accession_number}/{accession_number}_2.fastq")

def trim_reads(accession_number, read_type, fastp_path):
    """Trim reads using fastp."""
    # Skip trimming if already trimmed files exist
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq")):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping trimming.", level="info")
        return

    if read_type == '1':
        log_message(f"\nTrimming {accession_number} with fastp (single-end)...", level="info")
        trim_command = f"{fastp_path} -i {accession_number}/{accession_number}.fastq -o {accession_number}/{accession_number}_trimmed.fq --thread=4"
        run_command(trim_command)
    elif read_type == '2':
        log_message(f"\nTrimming {accession_number} with fastp (paired-end)...", level="info")
        trim_command = f"{fastp_path} -i {accession_number}/{accession_number}_1.fastq -I {accession_number}/{accession_number}_2.fastq -o {accession_number}/{accession_number}_1_trimmed.fq -O {accession_number}/{accession_number}_2_trimmed.fq --thread=4"
        run_command(trim_command)

    shutil.move("fastp.html", f"{accession_number}/fastp.html")
    shutil.move("fastp.json", f"{accession_number}/fastp.json")
