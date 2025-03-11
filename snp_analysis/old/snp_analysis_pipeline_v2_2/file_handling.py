import os
from typing import List
from snp_analysis_pipeline_v2_2.logging_module import log_message

def read_accession_numbers(file_path: str) -> List[str]:
    """
    Read accession numbers from a provided text file.

    Args:
        file_path (str): Path to the accession list file.

    Returns:
        List[str]: List of accession numbers.

    Raises:
        FileNotFoundError: If the file is not found.
    """
    try:
        with open(file_path, "r") as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
        return accession_numbers
    except FileNotFoundError:
        log_message(f"The file {file_path} was not found.", level="error")
        raise

def is_file_empty(file_path: str) -> bool:
    """
    Check if a file exists and is empty.

    Args:
        file_path (str): Path to the file.

    Returns:
        bool: True if the file exists and is empty, False otherwise.
    """
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0

def delete_intermediate_files(accession_number: str, chromosome: str) -> None:
    """
    Delete intermediate files such as SAM, BAM, and raw BCF files.

    Args:
        accession_number (str): The accession number of the sample.
        chromosome (str): The chromosome being processed.
    """
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

def detect_accession_list_file() -> str:
    """
    Detects a .txt file in the current directory and asks the user to select the correct one.

    Returns:
        str: Path to the selected accession list file.

    Raises:
        FileNotFoundError: If no valid file is selected.
    """
    txt_files = [f for f in os.listdir() if f.endswith(".txt")]

    if not txt_files:
        log_message("No .txt files found in the current directory. Please create an accession list file.", level="error")
        raise FileNotFoundError("No .txt files found. Please create an accession list file (one accession per line).")

    for txt_file in txt_files:
        log_message(f"Detected accession list file: {txt_file}. Would you like to use this file? (yes/y or no/n)", level="info")
        user_input = input().strip().lower()

        if user_input in ["yes", "y"]:
            return txt_file
        else:
            log_message(f"{txt_file} is not the correct file. Checking for another file...", level="warning")

    # If the loop completes without finding a valid file
    log_message("No valid accession list file was selected. Please ensure you have a correct file.", level="error")
    raise FileNotFoundError("No valid accession list file selected.")
