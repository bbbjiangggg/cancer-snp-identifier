import os
import subprocess
from pathlib import Path

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def run_analysis(sra):
    cwd = Path.cwd()
    print(f'{MAGENTA}Running SNP Analysis for {sra}...{RESET}')

    # Add your analysis logic here

    print(f'{GREEN}Analysis for {sra} Completed!{RESET}')

if __name__ == "__main__":
    # Example usage with a command-line argument
    import sys
    if len(sys.argv) > 1:
        sra = sys.argv[1]
        run_analysis(sra)
    else:
        print(f'{RED}Error: Please provide an SRA argument.{RESET}')
