import os
from snp_analysis.snp_analysis_pipeline_v2.logging_module import log_message

def read_accession_numbers(file_path):
    """Read accession numbers from a provided text file."""
    try:
        with open(file_path, 'r') as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        log_message(f"The file {file_path} was not found.", level="error")
        raise

def is_file_empty(file_path):
    """Check if a file exists and is empty."""
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0

def delete_intermediate_files(accession_number, chromosome):
    """Delete intermediate files such as SAM, BAM, and raw BCF files."""
    intermediate_files = [
        f"{accession_number}/{accession_number}_mapped_{chromosome}.sam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.bam",
        f"{accession_number}/{accession_number}_mapped_{chromosome}.raw.bcf"
    ]
    
    for file_path in intermediate_files:
        if os.path.isfile(file_path):
            os.remove(file_path)
            log_message(f"Deleted {file_path}", level="success")
        else:
            log_message(f"File {file_path} not found. Skipping deletion.", level="warning")

def detect_accession_list_file():
    """Detects a .txt file in the current directory and asks the user if it's the correct one."""
    txt_files = [f for f in os.listdir() if f.endswith('.txt')]
    if not txt_files:
        raise FileNotFoundError("No .txt files found in the current directory.")
    
    for txt_file in txt_files:
        log_message(f"Detected accession list file: {txt_file}. Is this correct? (yes/y or no/n)", level="info")
        user_input = input().strip().lower()
        if user_input in ['yes', 'y']:
            return txt_file
        else:
            log_message(f"{txt_file} is not the correct file. Moving to the next one if available.", level="warning")

    # If no valid file was selected, raise an error
    raise FileNotFoundError("No valid accession list file selected.")
