# path_management.py
import os
from snp_analysis.snp_analysis_pipeline_v2.logging_module import log_message

def ensure_directory(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        log_message(f"Directory {path} created.", level="success")
    else:
        log_message(f"Directory {path} already exists.", level="warning")

def print_chromosome_paths(chromosomes_list, bwa_base_path, bowtie_base_path):
    for chromosome in chromosomes_list:
        bwa_chrom_path = f"{bwa_base_path}hg38.fa"
        bowtie_index_path = f"{bowtie_base_path}hg38_index"
        log_message(f"Paths for chromosome {chromosome}:", level="info")
        log_message("BWA Chromosome Path: ", normal_part=bwa_chrom_path, level="info")
        log_message("Bowtie Index Path: ", normal_part=bowtie_index_path, level="info")
