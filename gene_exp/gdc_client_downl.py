import os
import subprocess

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Define function to run commands
def run_command(cmd):
    print(f'{MAGENTA}Running command: {cmd}{RESET}')
    subprocess.run(cmd, shell=True, check=True)

# Define paths and parameters
gdc_client_path = subprocess.check_output(['which', 'gdc-client']).strip().decode()
reference_genome_path = input(f'{MAGENTA}Enter the path to the reference genome in FASTA format:{RESET} ')
hisat2_index_path = input(f'{MAGENTA}Enter the path to the HISAT2 index:{RESET} ')
output_dir = input(f'{MAGENTA}Enter the path to the output directory:{RESET} ')

# Check if HISAT2 is installed
try:
    run_command('hisat2 --version')
except FileNotFoundError:
    print(f'{RED}HISAT2 is not installed. Please install it first.{RESET}')
    exit()

# Check if GDC Client is installed
try:
    run_command(f'{gdc_client_path} version')
except FileNotFoundError:
    print(f'{RED}GDC Client is not installed. Please install it first.{RESET}')
    exit()

# Download data from GDC
data_list_file = input(f'{MAGENTA}Enter the path to the data list file:{RESET} ')
with open(data_list_file, 'r') as f:
    for line in f:
        line = line.strip()
        data_id, bam_file, _, _, _, _, _ = line.split('\t')
        output_file_name = os.path.splitext(bam_file)[0]
        output_dir_name = os.path.join(output_dir, output_file_name)
        os.makedirs(output_dir_name, exist_ok=True)
        output_bam_path = os.path.join(output_dir_name, output_file_name + '_sorted.bam')
        
        # Check if BAM file already exists
        if os.path.exists(output_bam_path):
            print(f'{BLUE}{output_bam_path} already exists. Skipping...{RESET}')
            continue
        
        download_cmd = f'{gdc_client_path} download -m {data_id} -t RNA-Seq -d {output_dir_name}'
        run_command(download_cmd)

        # Map reads to reference genome using HISAT2
        hisat2_cmd = f'hisat2 -p 4 -x {hisat2_index_path} -U {os.path.join(output_dir_name, bam_file)} | samtools view -bS - | samtools sort -o {output_bam_path} -'
        run_command(hisat2_cmd)

        # Index the BAM file
        index_cmd = f'samtools index {output_bam_path}'
        run_command(index_cmd)

        # Run samtools flagstat to check the alignment
        flagstat = subprocess.check_output(['samtools', 'flagstat', output_bam_path])
        if b'error' in flagstat:
            print(f'{RED}Error in alignment of {output_bam_path}{RESET}')
            continue
        else:
            print(f'{GREEN}{output_bam_path} successfully generated!{RESET}')
