#!/usr/bin/env python3

import os

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'


# Set the parent directory path
parent_directory_path = os.getcwd()
print(f"{MAGENTA}\nParent directory path:{RESET} {parent_directory_path} ")

print(f"{MAGENTA}\nChecking subdirectories for var.-final.vcf...{RESET}")

# Get a list of all the subdirectories in the parent directory path
subdirectories = [d for d in os.listdir(parent_directory_path) if os.path.isdir(os.path.join(parent_directory_path, d))]

# Loop through each subdirectory
for subdir in subdirectories:
    # Check if the subdirectory starts with SRR or ERR
    if subdir.startswith("SRR") or subdir.startswith("ERR"):
        # Get the full path to the subdirectory
        subdirectory_path = os.path.join(parent_directory_path, subdir)
        # Get a list of all the files in the subdirectory
        files = os.listdir(subdirectory_path)
        # Check if any of the files end with "var.-final.vcf"
        if any(file.endswith("var.-final.vcf") for file in files):
            # If so, continue
            continue
        else:
            # If not, print a message saying so
            subdirectory_name = os.path.basename(subdirectory_path)
            print(f"{MAGENTA}\nNo file ending with 'var.-final.vcf' found in subdirectory:{RESET} {subdirectory_name}")
            print(f"{MAGENTA}\nIt needs to be removed or reanalyzed.{RESET}")
            # Give option to remove
            remove = input(f"{BLUE}\nRemove? (yes/no){RESET} ")
            if remove.lower() == "yes" or remove.lower() == "y":
                # Remove the subdirectory and its contents
                os.system(f"rm -r {subdirectory_path}")
                print(f"{MAGENTA}\nRemoved subdirectory{RESET} {subdirectory_path}.")
            else:
                print(f"{MAGENTA}\nSubdirectory not removed.{RESET}")





