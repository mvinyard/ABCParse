# -- import local modules: ----------------------------------------------------
from ._abc_logger import ABCLogger

# -- global logger instance: --------------------------------------------------
_default_logger = ABCLogger()

# -- import API-facing functions: ---------------------------------------------
from ._get_logger import get_logger
from ._set_global_log_level import set_global_log_level
from ._convenience_functions import debug, info, warning, error, critical
from . import _format

__all__ = [
    "ABCLogger",
    "_default_logger",
    "get_logger",
    "set_global_log_level",
    "debug",
    "info",
    "warning",
    "_format",
]