# -- import packages: ---------------------------------------------------------
import inspect

# -- set type hints: ----------------------------------------------------------
from typing import Any, Callable, Dict, List, Union

# -- import internal modules: -------------------------------------------------
from . import logging

# -- module logger: -----------------------------------------------------------
_logger = logging.get_logger(name="function_kwargs")


# -- operational class: -------------------------------------------------------
class KwargExtractor:
    def __init__(self, func):
        self.func = func
        self._logger = logging.get_logger(name="KwargExtractor", level="warning")
        self._logger.debug(f"Initializing KwargExtractor for function: {func.__name__}")

    # -- methods: -------------------------------------------------------------
    def __parse__(self, kwargs, parse_ignore=["self"]) -> None:
        self._PASSED = {}

        if isinstance(kwargs, dict):
            self._logger.debug(f"Parsing dictionary kwargs with {len(kwargs)} items")
            for key, val in kwargs.items():
                if not key in self._ignore:
                    self._PASSED[key] = val
                    self._logger.debug(f"Added key to _PASSED: {key}")
        elif isinstance(kwargs, list):
            self._logger.debug(f"Parsing list kwargs with {len(kwargs)} items")
            for key in kwargs:
                self._PASSED[key] = None
                self._logger.debug(f"Added key to _PASSED with None value: {key}")

    def _extract_func_params(self):
        params = list(inspect.signature(self.func).parameters.keys())
        self._logger.debug(f"Extracted function parameters: {params}")
        return params

    def query(self):
        self.func_kwargs = {}
        self._logger.debug("Querying for function kwargs")

        for key, val in self._PASSED.items():
            if not key in self._ignore:
                if key in self.func_params:
                    if self._obj:
                        val = getattr(self._obj, key)
                        self._logger.debug(f"Got attribute {key} from object")
                    self.func_kwargs[key] = val
                    self._logger.debug(f"Added to func_kwargs: {key}")

            if key == "kwargs":
                if not val is None:
                    self._check_literal_kwargs(val)

        self._logger.debug(f"Extracted {len(self.func_kwargs)} kwargs for function {self.func.__name__}")
        return self.func_kwargs

    def _check_literal_kwargs(self, kwargs):
        self._logger.debug(f"Checking literal kwargs: {kwargs}")
        for key, val in kwargs.items():
            if not key in self._ignore:
                self.func_kwargs[key] = val
                self._logger.debug(f"Added literal kwarg to func_kwargs: {key}")

    # -- property: ------------------------------------------------------------
    @property
    def func_params(self):
        return self._extract_func_params()

    def __call__(self, kwargs: Union[Dict, List], obj=None, ignore=["self", "kwargs"]):
        self._obj = obj
        self._logger.debug(f"Called with kwargs: {kwargs}, obj: {obj}, ignore: {ignore}")

        if (
            (not self._obj is None)
            and (not isinstance(kwargs, dict))
            and (not isinstance(kwargs, list))
        ):
            self._logger.debug(f"Using object directory as kwargs")
            kwargs = self._obj.__dir__()

        self._ignore = ignore
        self.__parse__(kwargs, parse_ignore=["self"])

        return self.query()


# -- API-facing function: -----------------------------------------------------
def function_kwargs(
    func: Callable,
    kwargs: Dict[str,Any] = None,
    obj: Any = None,
    ignore: List[str] = ["self", "kwargs"],
):
    """
    Returns the subset of kwargs that can be used in the func.

    Args:
        func (Callable): 
    
        kwargs
            if obj is passed, this argument is overridden.
    
        obj
            if kwargs is passed, obj overrides.

    Returns:
        function_kwargs
            type: list
    """
    _logger.debug(f"function_kwargs called for function: {func.__name__}")
    kwarg_extractor = KwargExtractor(func=func)
    result = kwarg_extractor(kwargs=kwargs, obj=obj, ignore=ignore)
    _logger.debug(f"Extracted {len(result)} kwargs for function {func.__name__}")
    return result
