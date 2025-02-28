# ABCParse

![Python Tests](https://github.com/mvinyard/ABCParse/actions/workflows/python-tests.yml/badge.svg)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/ABCParse.svg)](https://pypi.python.org/pypi/ABCParse/)
[![PyPI version](https://badge.fury.io/py/ABCParse.svg)](https://badge.fury.io/py/ABCParse)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A better base class that handles parsing local arguments.

```bash
pip install ABCParse
```

```python
from ABCParse import ABCParse


class SomeClass(ABCParse):
    def __init__(self, arg1, arg2):
      self.__parse__(kwargs=locals())
      
something = SomeClass(arg1 = 4, arg2 = "name")
```

## Logging

ABCParse includes a built-in logging system that provides visibility into the parsing process. You can configure the logging level, format, and output destination.

### Basic Usage

```python
from ABCParse import set_global_log_level, get_logger

# Set the global log level
set_global_log_level("debug")  # Options: debug, info, warning, error, critical

# Get a custom logger for your module
logger = get_logger(name="my_module")
logger.info("This is an info message")
logger.debug("This is a debug message")
```

### Advanced Configuration

```python
from ABCParse import get_logger

# Configure a logger with custom settings
logger = get_logger(
    name="my_custom_logger",
    level="debug",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    date_format="%Y-%m-%d %H:%M:%S",
    file_path="logs/my_app.log",  # Optional: log to file
    propagate=False
)

# Log messages at different levels
logger.debug("Detailed information for debugging")
logger.info("General information about program execution")
logger.warning("Warning about potential issues")
logger.error("Error that occurred during execution")
logger.critical("Critical error that may cause program failure")

# Log dictionary data
config = {"param1": 42, "param2": "value", "enabled": True}
logger.log_dict(config, level="debug", prefix="Configuration")
```
