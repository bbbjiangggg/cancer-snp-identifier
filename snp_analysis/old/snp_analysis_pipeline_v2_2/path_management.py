# path_management.py
import os
from typing import List
from snp_analysis_pipeline_v2_2.logging_module import log_message

def ensure_directory(path: str) -> None:
    """
    Ensure a directory exists. If not, create it.

    Args:
        path (str): Path to the directory.
    """
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        log_message(f"Directory {path} created.", level="success")
    else:
        log_message(f"Directory {path} already exists.", level="warning")

def print_chromosome_paths(chromosomes_list: List[str], bwa_base_path: str, bowtie_base_path: str) -> None:
    """
    Print paths for BWA and Bowtie2 indices for each chromosome.

    Args:
        chromosomes_list (List[str]): List of chromosomes to analyze.
        bwa_base_path (str): Base path for BWA indices.
        bowtie_base_path (str): Base path for Bowtie2 indices.
    """
    for chromosome in chromosomes_list:
        bwa_chrom_path = f"{bwa_base_path}hg38.fa"
        bowtie_index_path = f"{bowtie_base_path}hg38_index"
        log_message(f"Paths for chromosome {chromosome}:", level="info")
        log_message("BWA Chromosome Path: ", normal_part=bwa_chrom_path, level="info")
        log_message("Bowtie Index Path: ", normal_part=bowtie_index_path, level="info")