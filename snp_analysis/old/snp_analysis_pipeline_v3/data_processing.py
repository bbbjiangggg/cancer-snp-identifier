import os
import shutil
from snp_analysis_pipeline_v3.logging_module import log_message
from snp_analysis_pipeline_v3.command_execution import run_command

def detect_read_type(accession_number):
    """
    Automatically detect if reads are single-end or paired-end.
    This function is no longer needed for the simplified script, but it can be kept for future use.
    """
    pass

def prefetch_and_convert(accession_number, threads):
    """
    Download and convert SRA files using fastq-dump without compression.
    This function only downloads the SRR file and moves it to the appropriate directory.
    """
    # Ensure the directory for the accession number exists
    os.makedirs(accession_number, exist_ok=True)

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
    """
    This function is no longer needed for the simplified script.
    """
    pass