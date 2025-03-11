# data_processing.py
import os
import shutil
from typing import Optional
from snp_analysis_pipeline_v2_2.logging_module import log_message
from snp_analysis_pipeline_v2_2.command_execution import run_command

def detect_read_type(accession_number: str) -> str:
    """
    Automatically detect if reads are single-end or paired-end.

    Args:
        accession_number (str): The accession number of the sample.

    Returns:
        str: '1' for single-end, '2' for paired-end.

    Raises:
        FileNotFoundError: If no valid FASTQ files are found.
    """
    # Check if trimmed files exist and skip FASTQ detection
    trimmed_single = f"{accession_number}/{accession_number}_trimmed.fq"
    trimmed_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
    trimmed_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

    if os.path.isfile(trimmed_single) or (os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2)):
        log_message(f"Trimmed files detected for {accession_number}. Skipping FASTQ detection.", level="info")
        if os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2):
            return '2'
        else:
            return '1'

    # If no trimmed files are found, check for FASTQ files
    fastq_single = f"{accession_number}/{accession_number}.fastq"
    fastq_paired_1 = f"{accession_number}/{accession_number}_1.fastq"
    fastq_paired_2 = f"{accession_number}/{accession_number}_2.fastq"

    if os.path.isfile(fastq_paired_1) and os.path.isfile(fastq_paired_2):
        log_message(f"Detected paired-end reads for {accession_number}.", level="info")
        return '2'
    elif os.path.isfile(fastq_single):
        log_message(f"Detected single-end reads for {accession_number}.", level="info")
        return '1'
    else:
        log_message(f"No valid FASTQ files found for {accession_number}.", level="error")
        raise FileNotFoundError(f"No FASTQ files found for {accession_number}.")

def prefetch_and_convert(accession_number: str, threads: int, prefetch_path: str, fasterq_dump_path: str) -> None:
    """
    Run prefetch and fasterq-dump to download and convert SRA files.

    Args:
        accession_number (str): The accession number of the sample.
        threads (int): Number of threads to use for fasterq-dump.
        prefetch_path (str): Path to the prefetch executable.
        fasterq_dump_path (str): Path to the fasterq-dump executable.

    Raises:
        FileNotFoundError: If the prefetch or fasterq-dump commands fail.
    """
    # Skip downloading if trimmed files are already present
    trimmed_single = f"{accession_number}/{accession_number}_trimmed.fq"
    trimmed_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
    trimmed_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

    if os.path.isfile(trimmed_single) or (os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2)):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping download and conversion.", level="info")
        return

    # Download the SRA file
    log_message(f"Downloading sequence number {accession_number} with prefetch...", level="info")
    run_command(f"{prefetch_path} {accession_number}")

    # Convert the SRA file to FASTQ
    log_message(f"Converting {accession_number} with fasterq-dump...", level="info")
    run_command(f"{fasterq_dump_path} {accession_number} --threads {threads} --progress")

    # Move the FASTQ files to the accession directory
    if os.path.isfile(f"{accession_number}.fastq"):
        shutil.move(f"{accession_number}.fastq", f"{accession_number}/{accession_number}.fastq")
    if os.path.isfile(f"{accession_number}_1.fastq"):
        shutil.move