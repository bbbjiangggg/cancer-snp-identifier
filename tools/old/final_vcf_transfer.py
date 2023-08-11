import os
import shutil

# Prompt the user to enter the path to the directory
dir_path = input("Enter the path to the directory to search for files in: ")

# Loop through all files in the directory
for filename in os.listdir(dir_path):
    # Check if the file name starts with "SRR" or "ERR"
    if filename.startswith(('SRR', 'ERR')):
        # Get the matching subdirectory name
        dir_name = filename.split('_')[0]

        # Check if the subdirectory exists in the current directory
        if os.path.isdir(os.path.join(dir_path, dir_name)):
            # Move the file to the matching subdirectory
            if filename.endswith('final.vcf'):
                shutil.move(os.path.join(dir_path, filename), os.path.join(dir_path, dir_name, filename))
        else:
            print(f"The directory {dir_name} does not exist in {dir_path}.")

