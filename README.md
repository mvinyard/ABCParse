# ABCParse

![Python Tests](https://github.com/mvinyard/AutoParser/actions/workflows/python-tests.yml/badge.svg)

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
