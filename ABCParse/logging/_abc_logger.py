# -- import packages: ---------------------------------------------------------
import logging
import pathlib
import sys

# -- import local modules: ----------------------------------------------------
from ._format import DEFAULT_LOG_FORMAT, DEFAULT_DATE_FORMAT, DEFAULT_LOG_FILEPATH

# -- set type hints: ----------------------------------------------------------
from typing import Any, Dict, Optional

# -- set constants: -----------------------------------------------------------
LOG_LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL
}

# -- operational class: -------------------------------------------------------
class ABCLogger:
    """
    Configurable logger for ABCParse.
    
    This class provides a centralized logging system for the ABCParse package,
    allowing for consistent logging across all modules.
    """
    
    def __init__(
        self,
        name: str = "ABCParse",
        level: str = "warning",
        format: str = DEFAULT_LOG_FORMAT,
        date_format: str = DEFAULT_DATE_FORMAT,
        file_path: Optional[str] = DEFAULT_LOG_FILEPATH,
        propagate: bool = False,
        file_level: Optional[str] = "debug",
    ) -> None:
        """
        Initialize the logger.
        
        Parameters
        ----------
        name : str, default="ABCParse"
            Name of the logger.
        level : str, default="info"
            Logging level for the console handler. One of: "debug", "info", "warning", "error", "critical".
        format : str, default=DEFAULT_LOG_FORMAT
            Log message format.
        date_format : str, default=DEFAULT_DATE_FORMAT
            Date format for log messages.
        file_path : Optional[str], default=None
            If provided, logs will be written to this file in addition to stdout.
        propagate : bool, default=False
            Whether to propagate logs to parent loggers.
        file_level : Optional[str], default=None
            Logging level for the file handler. If None, uses the same level as the console handler.
        """

        self.name = name
        self.level = level
        self.file_level = file_level if file_level else level
        self.format = format
        self.date_format = date_format
        self.file_path = file_path
        self.file_handler = None
        self._file = None
        self.console_handler = None
        
        # Create logger
        self.logger = logging.getLogger(name)
        self.logger.setLevel(min(
            LOG_LEVELS.get(level.lower(), logging.INFO),
            LOG_LEVELS.get(self.file_level.lower(), logging.INFO)
        ))
        self.logger.propagate = propagate
        
        # Remove existing handlers to avoid duplicates when reconfiguring
        self._remove_handlers()
        
        # Create formatter
        self.formatter = logging.Formatter(format, date_format)
        
        # Create console handler
        self.console_handler = logging.StreamHandler(sys.stdout)
        self.console_handler.setFormatter(self.formatter)
        self.console_handler.setLevel(LOG_LEVELS.get(level.lower(), logging.INFO))
        self.logger.addHandler(self.console_handler)
        
        # Create file handler if file_path is provided
        if file_path:
            self._setup_file_logging(file_path)
    
    def _remove_handlers(self) -> None:
        """Remove and close all existing handlers."""
        if hasattr(self, 'logger') and self.logger and self.logger.handlers:
            for handler in list(self.logger.handlers):
                try:
                    if isinstance(handler, logging.FileHandler):
                        handler.close()
                    self.logger.removeHandler(handler)
                except Exception as e:
                    print(f"Error removing handler: {e}")
    
    def _setup_file_logging(self, file_path) -> None:
        """Set up file logging with proper error handling."""
        try:
            # Ensure directory exists
            if not pathlib.Path(file_path).parent.exists():
                pathlib.Path(file_path).parent.mkdir(parents=True, exist_ok=True)
                        
            # Open the file directly for more control
            self._file = open(file_path, 'w')
            
            # Also create a standard file handler
            self.file_handler = logging.FileHandler(file_path)
            self.file_handler.setFormatter(self.formatter)
            self.file_handler.setLevel(LOG_LEVELS.get(self.file_level.lower(), logging.INFO))
            self.logger.addHandler(self.file_handler)
        except Exception as e:
            print(f"Error setting up file logging: {e}")
            # Clean up any partially initialized resources
            self._close_file_resources()
            
    def _write_to_file(self, msg) -> None:
        """Write directly to the file if it exists."""
        if self._file and not self._file.closed:
            try:
                self._file.write(msg + '\n')
                self._file.flush()
            except Exception as e:
                print(f"Error writing to file: {e}")
                self._close_file_resources()
    
    def _close_file_resources(self) -> None:
        """Close file resources with proper error handling."""
        # Close the direct file if it exists
        if self._file:
            try:
                if not self._file.closed:
                    self._file.flush()
                    self._file.close()
            except Exception as e:
                print(f"Error closing file: {e}")
            finally:
                self._file = None
        
        # Close and remove the file handler if it exists
        if self.file_handler:
            try:
                self.file_handler.close()
                if self.logger:
                    self.logger.removeHandler(self.file_handler)
            except Exception as e:
                print(f"Error closing file handler: {e}")
            finally:
                self.file_handler = None
            
    def close(self) -> None:
        """
        Close all handlers associated with this logger.
        This is particularly important for file handlers to ensure they flush and close properly.
        """
        if hasattr(self, 'logger') and self.logger:
            for handler in list(self.logger.handlers):
                try:
                    handler.flush()
                    if isinstance(handler, logging.FileHandler):
                        handler.close()
                        self.logger.removeHandler(handler)
                except Exception as e:
                    print(f"Error closing handler: {e}")
        
        # Close file resources
        self._close_file_resources()
                
    def __del__(self) -> None:
        """Ensure file handlers are closed when the logger is garbage collected."""
        try:
            self.close()
        except:
            pass
    
    def set_level(self, level: str, handler_type: str = "both") -> None:
        """
        Set the logging level.
        
        Parameters
        ----------
        level : str
            Logging level. One of: "debug", "info", "warning", "error", "critical".
        handler_type : str, default="both"
            Which handler to set the level for. One of: "console", "file", "both".
        """
        if level.lower() in LOG_LEVELS:
            log_level = LOG_LEVELS[level.lower()]
            
            if handler_type.lower() in ["console", "both"]:
                if self.console_handler:
                    self.console_handler.setLevel(log_level)
                self.level = level.lower()
                
            if handler_type.lower() in ["file", "both"]:
                if self.file_handler:
                    self.file_handler.setLevel(log_level)
                self.file_level = level.lower()
                
            # Set the logger's level to the minimum of both handlers
            min_level = min(
                LOG_LEVELS.get(self.level.lower(), logging.INFO),
                LOG_LEVELS.get(self.file_level.lower(), logging.INFO)
            )
            self.logger.setLevel(min_level)
    
    def debug(self, msg: str, *args, **kwargs) -> None:
        """Log a debug message."""
        self.logger.debug(msg, *args, **kwargs)
        if LOG_LEVELS.get(self.file_level.lower(), logging.INFO) <= logging.DEBUG:
            self._write_to_file(f"DEBUG: {msg}")
        sys.stdout.flush()
    
    def info(self, msg: str, *args, **kwargs) -> None:
        """Log an info message."""
        self.logger.info(msg, *args, **kwargs)
        if LOG_LEVELS.get(self.file_level.lower(), logging.INFO) <= logging.INFO:
            self._write_to_file(f"INFO: {msg}")
        sys.stdout.flush()
    
    def warning(self, msg: str, *args, **kwargs) -> None:
        """Log a warning message."""
        self.logger.warning(msg, *args, **kwargs)
        if LOG_LEVELS.get(self.file_level.lower(), logging.INFO) <= logging.WARNING:
            self._write_to_file(f"WARNING: {msg}")
        sys.stdout.flush()
    
    def error(self, msg: str, *args, **kwargs) -> None:
        """Log an error message."""
        self.logger.error(msg, *args, **kwargs)
        if LOG_LEVELS.get(self.file_level.lower(), logging.INFO) <= logging.ERROR:
            self._write_to_file(f"ERROR: {msg}")
        sys.stdout.flush()
    
    def critical(self, msg: str, *args, **kwargs) -> None:
        """Log a critical message."""
        self.logger.critical(msg, *args, **kwargs)
        if LOG_LEVELS.get(self.file_level.lower(), logging.INFO) <= logging.CRITICAL:
            self._write_to_file(f"CRITICAL: {msg}")
        sys.stdout.flush()
    
    def log_dict(self, data: Dict[str, Any], level: str = "info", prefix: str = "") -> None:
        """
        Log a dictionary with key-value pairs.
        
        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary to log.
        level : str, default="info"
            Logging level for the messages.
        prefix : str, default=""
            Prefix to add to each log message.
        """
        log_method = getattr(self, level.lower(), self.info)
        
        for key, value in data.items():
            if prefix:
                log_method(f"{prefix} - {key}: {value}")
            else:
                log_method(f"{key}: {value}")