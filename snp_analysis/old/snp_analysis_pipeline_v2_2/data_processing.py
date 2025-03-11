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
    trimmed_single = f"{accession_number}/{accession_number}_trimmed.fq"
    trimmed_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
    trimmed_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

    if os.path.isfile(trimmed_single) or (os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2)):
        log_message(f"Trimmed files detected for {accession_number}. Skipping FASTQ detection.", level="info")
        return '2' if os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2) else '1'

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
    """
    trimmed_single = f"{accession_number}/{accession_number}_trimmed.fq"
    trimmed_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
    trimmed_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

    if os.path.isfile(trimmed_single) or (os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2)):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping download and conversion.", level="info")
        return

    log_message(f"Downloading sequence number {accession_number} with prefetch...", level="info")
    run_command(f"{prefetch_path} {accession_number}")

    log_message(f"Converting {accession_number} with fasterq-dump...", level="info")
    run_command(f"{fasterq_dump_path} {accession_number} --threads {threads} --progress")

    if os.path.isfile(f"{accession_number}.fastq"):
        shutil.move(f"{accession_number}.fastq", f"{accession_number}/{accession_number}.fastq")
    if os.path.isfile(f"{accession_number}_1.fastq"):
        shutil.move(f"{accession_number}_1.fastq", f"{accession_number}/{accession_number}_1.fastq")
    if os.path.isfile(f"{accession_number}_2.fastq"):
        shutil.move(f"{accession_number}_2.fastq", f"{accession_number}/{accession_number}_2.fastq")

def trim_reads(accession_number: str, read_type: str, fastp_path: str) -> None:
    """
    Trim reads using fastp.
    """
    trimmed_single = f"{accession_number}/{accession_number}_trimmed.fq"
    trimmed_paired_1 = f"{accession_number}/{accession_number}_1_trimmed.fq"
    trimmed_paired_2 = f"{accession_number}/{accession_number}_2_trimmed.fq"

    if os.path.isfile(trimmed_single) or (os.path.isfile(trimmed_paired_1) and os.path.isfile(trimmed_paired_2)):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping trimming.", level="info")
        return

    if read_type == '1':
        log_message(f"Trimming {accession_number} with fastp (single-end)...", level="info")
        trim_command = f"{fastp_path} -i {accession_number}/{accession_number}.fastq -o {accession_number}/{accession_number}_trimmed.fq --thread=4"
        run_command(trim_command)
    elif read_type == '2':
        log_message(f"Trimming {accession_number} with fastp (paired-end)...", level="info")
        trim_command = f"{fastp_path} -i {accession_number}/{accession_number}_1.fastq -I {accession_number}/{accession_number}_2.fastq -o {accession_number}/{accession_number}_1_trimmed.fq -O {accession_number}/{accession_number}_2_trimmed.fq --thread=4"
        run_command(trim_command)

    shutil.move("fastp.html", f"{accession_number}/fastp.html")
    shutil.move("fastp.json", f"{accession_number}/fastp.json")
