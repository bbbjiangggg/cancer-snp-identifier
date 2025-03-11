import os
import shutil
from snp_analysis_pipeline_v2_1.logging_module import log_message
from snp_analysis_pipeline_v2_1.command_execution import run_command

def detect_read_type(accession_number):
    """Automatically detect if reads are single-end or paired-end."""
    # Check if trimmed files exist and skip FASTQ detection
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq.gz") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq.gz") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq.gz")):
        log_message(f"Trimmed files detected for {accession_number}. Skipping FASTQ detection.", level="info")
        if os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq.gz") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq.gz"):
            return '2'
        else:
            return '1'

    # Check for compressed FASTQ files first
    if os.path.isfile(f"{accession_number}/{accession_number}_1.fastq.gz") and os.path.isfile(f"{accession_number}/{accession_number}_2.fastq.gz"):
        log_message(f"Detected paired-end reads for {accession_number}.", level="info")
        return '2'
    elif os.path.isfile(f"{accession_number}/{accession_number}.fastq.gz") or os.path.isfile(f"{accession_number}/{accession_number}_1.fastq.gz"):
        log_message(f"Detected single-end read for {accession_number}.", level="info")
        return '1'
    else:
        log_message(f"No valid FASTQ files found for {accession_number}.", level="error")
        raise FileNotFoundError(f"No FASTQ files found for {accession_number}.")

def prefetch_and_convert(accession_number, threads):
    """Download and convert SRA files using fastq-dump instead of prefetch + fasterq-dump."""
    # Skip downloading if trimmed files are already present
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq.gz") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq.gz") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq.gz")):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping download and conversion.", level="info")
        return

    log_message(f"\nDownloading and converting {accession_number} using fastq-dump...", level="info")

    # Run fastq-dump to download reads
    fastq_dump_command = f"fastq-dump --split-files --gzip {accession_number}"
    run_command(fastq_dump_command)

    # Check if output files exist and rename them to match the expected format
    single_end_file = f"{accession_number}.fastq.gz"
    paired_end_file_1 = f"{accession_number}_1.fastq.gz"
    paired_end_file_2 = f"{accession_number}_2.fastq.gz"

    # If paired-end files exist, keep them as they are
    if os.path.isfile(paired_end_file_1) and os.path.isfile(paired_end_file_2):
        log_message(f"Paired-end reads detected for {accession_number}.", level="info")
        shutil.move(paired_end_file_1, f"{accession_number}/{accession_number}_1.fastq.gz")
        shutil.move(paired_end_file_2, f"{accession_number}/{accession_number}_2.fastq.gz")
    
    # If only _1.fastq.gz exists, assume it's a single-end read and rename it
    elif os.path.isfile(paired_end_file_1) and not os.path.isfile(paired_end_file_2):
        log_message(f"Single-end read detected for {accession_number}. Renaming file.", level="info")
        shutil.move(paired_end_file_1, f"{accession_number}/{accession_number}.fastq.gz")

    # If the expected files are missing, log an error
    else:
        log_message(f"❌ FASTQ file(s) missing for {accession_number}. Check fastq-dump output.", level="error")
        raise FileNotFoundError(f"FASTQ file(s) missing for {accession_number}.")


def trim_reads(accession_number, read_type, fastp_path):
    """Trim reads using fastp."""
    # Skip trimming if already trimmed files exist
    if os.path.isfile(f"{accession_number}/{accession_number}_trimmed.fq.gz") or (
        os.path.isfile(f"{accession_number}/{accession_number}_1_trimmed.fq.gz") and os.path.isfile(f"{accession_number}/{accession_number}_2_trimmed.fq.gz")):
        log_message(f"Trimmed files already exist for {accession_number}. Skipping trimming.", level="info")
        return

    if read_type == '1':
        log_message(f"\nTrimming {accession_number} with fastp (single-end)...", level="info")
        trim_command = f"{fastp_path} -i {accession_number}/{accession_number}.fastq.gz -o {accession_number}/{accession_number}_trimmed.fq.gz --thread=4"
        run_command(trim_command)
    elif read_type == '2':
        log_message(f"\nTrimming {accession_number} with fastp (paired-end)...", level="info")
        trim_command = f"{fastp_path} -i {accession_number}/{accession_number}_1.fastq.gz -I {accession_number}/{accession_number}_2.fastq.gz -o {accession_number}/{accession_number}_1_trimmed.fq.gz -O {accession_number}/{accession_number}_2_trimmed.fq.gz --thread=4"
        run_command(trim_command)

    shutil.move("fastp.html", f"{accession_number}/fastp.html")
    shutil.move("fastp.json", f"{accession_number}/fastp.json")
