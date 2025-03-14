# command_execution.py
import subprocess
import sys
from snp_analysis_pipeline_v2_1.logging_module import log_message

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        log_message(f"An error occurred: {e}", level="error")
        sys.exit(1)
    except KeyboardInterrupt:
        log_message("Analysis interrupted by user. Exiting.", level="error")
        sys.exit(1)
