# -- import local modules: ----------------------------------------------------
from ._abc_logger import ABCLogger
from ._format import DEFAULT_LOG_FORMAT, DEFAULT_DATE_FORMAT

# -- set type hints: ----------------------------------------------------------
from typing import Optional

# -- API-facing functions: ----------------------------------------------------
def get_logger(
    name: str = "ABCParse",
    level: str = "info",
    format: str = DEFAULT_LOG_FORMAT,
    date_format: str = DEFAULT_DATE_FORMAT,
    file_path: Optional[str] = None,
    propagate: bool = False
) -> ABCLogger:
    """
    Get a configured logger instance.
    
    Parameters
    ----------
    name : str, default="ABCParse"
        Name of the logger.
    level : str, default="info"
        Logging level. One of: "debug", "info", "warning", "error", "critical".
    format : str, default=DEFAULT_LOG_FORMAT
        Log message format.
    date_format : str, default=DEFAULT_DATE_FORMAT
        Date format for log messages.
    file_path : Optional[str], default=None
        If provided, logs will be written to this file in addition to stdout.
    propagate : bool, default=False
        Whether to propagate logs to parent loggers.
        
    Returns
    -------
    ABCLogger
        Configured logger instance.
    """
    return ABCLogger(
        name=name,
        level=level,
        format=format,
        date_format=date_format,
        file_path=file_path,
        propagate=propagate
    )