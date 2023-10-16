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
