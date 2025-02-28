from . import _default_logger
__module_name__ = "_logging.py"
__doc__ = """Convenience functions for the default logger."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["mvinyard.ai@gmail.com"])


# -- functions: ---------------------------------------------------------------
def debug(msg: str, *args, **kwargs) -> None:
    """Log a debug message using the default logger."""
    _default_logger.debug(msg, *args, **kwargs)

def info(msg: str, *args, **kwargs) -> None:
    """Log an info message using the default logger."""
    _default_logger.info(msg, *args, **kwargs)

def warning(msg: str, *args, **kwargs) -> None:
    """Log a warning message using the default logger."""
    _default_logger.warning(msg, *args, **kwargs)

def error(msg: str, *args, **kwargs) -> None:
    """Log an error message using the default logger."""
    _default_logger.error(msg, *args, **kwargs)

def critical(msg: str, *args, **kwargs) -> None:
    """Log a critical message using the default logger."""
    _default_logger.critical(msg, *args, **kwargs) 