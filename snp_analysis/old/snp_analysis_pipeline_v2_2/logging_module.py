# logging_module.py
import logging
from logging.handlers import RotatingFileHandler
import os

def setup_logging(log_file: str = "pipeline.log", log_level: str = "INFO", max_log_size: int = 5242880, backup_count: int = 3) -> None:
    """
    Set up logging configuration.

    Args:
        log_file (str): Path to the log file. Default is "pipeline.log".
        log_level (str): Logging level (e.g., "INFO", "WARNING"). Default is "INFO".
        max_log_size (int): Maximum size of a log file in bytes before rotation. Default is 5MB.
        backup_count (int): Number of backup log files to keep. Default is 3.
    """
    # Ensure the log directory exists
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)

    # Configure logging
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),  # Print to console
            RotatingFileHandler(log_file, maxBytes=max_log_size, backupCount=backup_count)  # Rotating log file
        ]
    )

    # Add a custom log level for "SUCCESS"
    SUCCESS_LEVEL_NUM = 25
    logging.addLevelName(SUCCESS_LEVEL_NUM, "SUCCESS")

    def success(self, message, *args, **kwargs):
        if self.isEnabledFor(SUCCESS_LEVEL_NUM):
            self._log(SUCCESS_LEVEL_NUM, message, args, **kwargs)

    logging.Logger.success = success

def log_message(message: str, level: str = "info") -> None:
    """
    Logs messages with different severity levels.

    Args:
        message (str): The message to log.
        level (str): The log level ("info", "warning", "error", "success"). Default is "info".
    """
    levels = {
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "success": 25  # Custom SUCCESS level
    }

    logger = logging.getLogger(__name__)
    logger.log(levels.get(level, logging.INFO), message)