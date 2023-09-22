import os

def remove_files_from_directory(path):
    """
    Cycle through directories with filenames starting with 'SRR' 
    and remove files that end with .fastq and .fastqc.zip.
    """
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if filename.startswith("SRR"):
                if filename.endswith(".fastq") or filename.endswith(".fastqc.zip"):
                    file_path = os.path.join(dirpath, filename)
                    try:
                        os.remove(file_path)
                        print(f"Deleted: {file_path}")
                    except Exception as e:
                        print(f"Error deleting {file_path}. Reason: {e}")

if __name__ == "__main__":
    directory_path = input("Enter the path to the root directory: ")
    remove_files_from_directory(directory_path)

