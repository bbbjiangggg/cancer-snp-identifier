# command_execution.py
import subprocess
import sys
from typing import Optional
from snp_analysis_pipeline_v2_1.logging_module import log_message

def run_command(command: str, timeout: Optional[int] = None) -> None:
    """
    Execute a shell command with error handling and optional timeout.

    Args:
        command (str): The command to execute.
        timeout (Optional[int]): Maximum time (in seconds) to wait for the command to complete. Default is None.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
        subprocess.TimeoutExpired: If the command does not complete within the specified timeout.
        KeyboardInterrupt: If the user interrupts the command execution.
    """
    try:
        # Execute the command
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            timeout=timeout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Log the command output
        if result.stdout:
            log_message(f"Command output: {result.stdout}", level="info")
        if result.stderr:
            log_message(f"Command error output: {result.stderr}", level="warning")

    except subprocess.CalledProcessError as e:
        # Log the error and exit if the command fails
        log_message(f"Command failed with exit code {e.returncode}: {e.stderr}", level="error")
        sys.exit(1)
    except subprocess.TimeoutExpired as e:
        # Log the timeout error and exit
        log_message(f"Command timed out after {timeout} seconds: {e.stderr}", level="error")
        sys.exit(1)
    except KeyboardInterrupt:
        # Handle user interruption
        log_message("Command execution interrupted by user. Exiting.", level="error")
        sys.exit(1)

def run_command_with_progress(command: str, description: str = "Processing", duration_estimate: int = 100) -> None:
    """
    Execute a shell command with a progress bar.

    Args:
        command (str): The command to execute.
        description (str): Description for the progress bar. Default is "Processing".
        duration_estimate (int): Estimated duration of the command (in seconds). Default is 100.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
        KeyboardInterrupt: If the user interrupts the command execution.
    """
    from tqdm import tqdm
    import time

    try:
        # Start the command
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Display a progress bar
        with tqdm(total=duration_estimate, desc=description, unit="s", ncols=100) as pbar:
            while process.poll() is None:
                time.sleep(1)  # Simulate work being done
                pbar.update(1)

            # Wait for the command to complete
            stdout, stderr = process.communicate()

            # Log the command output
            if stdout:
                log_message(f"Command output: {stdout}", level="info")
            if stderr:
                log_message(f"Command error output: {stderr}", level="warning")

            # Check for errors
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, command, output=stdout, stderr=stderr)

            # Finish the progress bar
            pbar.update(duration_estimate - pbar.n)

    except subprocess.CalledProcessError as e:
        # Log the error and exit if the command fails
        log_message(f"Command failed with exit code {e.returncode}: {e.stderr}", level="error")
        sys.exit(1)
    except KeyboardInterrupt:
        # Handle user interruption
        log_message("Command execution interrupted by user. Exiting.", level="error")
        sys.exit(1)