import os
import shutil
import subprocess
from pathlib import Path

# Define the output directory
output_dir = "crispr_output"

# Check if output directory exists, create it if necessary
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Get the current working directory
cwd = Path.cwd()
print('\n')
print(f'{MAGENTA} This is your current working directory: {RESET}{cwd}\n')

# Define the directory to search in
search_dir = cwd

def find_unanalyzed_files(src_dir):
    unanalyzed_files = []
    for root, dirs, files in os.walk(src_dir):
        if 'RR' in root:
            for file in files:
                if file.endswith("_mapped.sorted.bam"):
                    src_file = os.path.join(root, file)
                    dst_file = os.path.join(output_dir, file)
                    if not os.path.exists(dst_file):
                        shutil.copy(src_file, dst_file)
                    break
            else:
                rr_dir = os.path.basename(root)
                unanalyzed_files.append(rr_dir)
    return sorted(list(set(unanalyzed_files)))


def print_unanalyzed_files(unanalyzed_files):
    print(f'{MAGENTA}\nThese are the unanalyzed files in the current directory:{RESET}')
    for f in unanalyzed_files:
        print(f'RR directory: {f}')
    print(f'{MAGENTA}There are{RESET} {len(unanalyzed_files)} {MAGENTA}RR directories with no _mapped.sorted.bam file.{RESET}')

# Find and print the sorted list of unanalyzed files
unanalyzed_files = find_unanalyzed_files(search_dir)
print_unanalyzed_files(unanalyzed_files)

# Ask to continue the analysis even if not all RR directories have the file
while True:
    choice = input(f'{MAGENTA}\nDo you want to continue the analysis? Enter yes or no: {RESET}')

    if choice.lower() == 'yes':
        # Analyze the files in the output directory with CRISPResso2
        target_output = open("target_output.txt", "w")
        for file in os.listdir(output_dir):
            if file.endswith("_mapped.sorted.bam"):
                file_path = os.path.join(output_dir, file)
                file_name = os.path.splitext(file)[0]
                command = f"CRISPResso -r1 -a GATCGGAAGAGCACACGTCT -g GACTGGTTCCAATTGAAAGTGTGACTGGTGCCAAGGAACTGCATTC -n {file_name} -f {file_path} -o {output_dir}/{file_name}"
                os.system(command)
                target_output.write(f"Results for {file_name}:\n")
                with open(f"{output_dir}/{file_name}/CRISPResso2_info.txt") as info_file:
                    info = info_file.read()
                    target_output.write(info)
        target_output.close()

        # Print the results
        print(f"{GREEN}\nThe analysis is complete. Check the target_output.txt file for the results.{RESET}")
        break
    elif choice.lower() == 'no':
        # Exit the program
        print(f'{MAGENTA}\nThe analysis was terminated. Goodbye.{RESET}')
