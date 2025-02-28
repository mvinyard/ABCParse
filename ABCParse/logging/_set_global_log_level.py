
from . import _default_logger
def set_global_log_level(level: str) -> None:
    """
    Set the log level for the default logger.
    
    Parameters
    ----------
    level : str
        Logging level. One of: "debug", "info", "warning", "error", "critical".
    """
    _default_logger.set_level(level)