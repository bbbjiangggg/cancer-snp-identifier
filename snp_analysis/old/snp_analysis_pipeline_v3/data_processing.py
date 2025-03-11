import os
import shutil
from snp_analysis_pipeline_v3.logging_module import log_message
from snp_analysis_pipeline_v3.command_execution import run_command

def detect_read_type(accession_number):
    """Automatically detect if reads are single-end or paired-end."""
    # Check if trimmed files exist and skip FASTQ detection
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq")):
        log_message(f"Trimmed files detected for {accession_number}. Skipping FASTQ detection.", level="info")
        if os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq"):
            return '2'  # Paired-end
        else:
            return '1'  # Single-end

    # Check for uncompressed FASTQ files
    if os.path.isfile(f"{accession_number}/{accession_number}_1.fastq") and os.path.isfile(f"{accession_number}/{accession_number}_2.fastq"):
        log_message(f"Detected paired-end reads for {accession_number}.", level="info")
        return '2'  # Paired-end
    elif os.path.isfile(f"{accession_number}/{accession_number}.fastq"):
        log_message(f"Detected single-end read for {accession_number}.", level="info")
        return '1'  # Single-end
    else:
        log_message(f"No valid FASTQ files found for {accession_number}.", level="error")
        raise FileNotFoundError(f"No FASTQ files found for {accession_number}.")

def prefetch_and_convert(accession_number, threads):
    """Download and convert SRA files using fastq-dump without compression."""
    # Skip downloading if trimmed files are already present
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq")):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping download and conversion.", level="info")
        return

    log_message(f"\nDownloading and converting {accession_number} using fastq-dump...", level="info")
    
    # Run fastq-dump without compression
    fastq_dump_command = f"fastq-dump --split-files {accession_number}"
    run_command(fastq_dump_command)

    # Move output FASTQ files to the designated directory
    if os.path.isfile(f"{accession_number}.fastq"):
        shutil.move(f"{accession_number}.fastq", f"{accession_number}/{accession_number}.fastq")
    if os.path.isfile(f"{accession_number}_1.fastq"):
        shutil.move(f"{accession_number}_1.fastq", f"{accession_number}/{accession_number}_1.fastq")
    if os.path.isfile(f"{accession_number}_2.fastq"):
        shutil.move(f"{accession_number}_2.fastq", f"{accession_number}/{accession_number}_2.fastq")

    log_message(f"Download and conversion completed for {accession_number}.", level="info")

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