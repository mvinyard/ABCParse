
# -- import packages: --------------------------------------------------------------------
from typing import Any, Callable, Dict, List, Union
import inspect


# -- operational class: ------------------------------------------------------------------
class KwargExtractor:
    def __init__(self, func):
        self.func = func

    # -- methods: ------------------------------------------------------------------------
    def __parse__(self, kwargs, parse_ignore=["self"]):
        self._PASSED = {}

        if isinstance(kwargs, dict):
            for key, val in kwargs.items():
                if not key in self._ignore:
                    self._PASSED[key] = val
        elif isinstance(kwargs, list):
            for key in kwargs:
                self._PASSED[key] = None

    def _extract_func_params(self):
        return list(inspect.signature(self.func).parameters.keys())

    def query(self):
        self.func_kwargs = {}

        for key, val in self._PASSED.items():
            if not key in self._ignore:
                if key in self.func_params:
                    if self._obj:
                        val = getattr(self._obj, key)
                    #                     print(key, "...regular")
                    self.func_kwargs[key] = val

            if key == "kwargs":
                if not val is None:
                    self._check_literal_kwargs(val)

        return self.func_kwargs

    def _check_literal_kwargs(self, kwargs):
        for key, val in kwargs.items():
            if not key in self._ignore:
                #                 print(key, "...literal kwargs")
                self.func_kwargs[key] = val

    # -- property: -----------------------------------------------------------------------
    @property
    def func_params(self):
        return self._extract_func_params()

    def __call__(self, kwargs: Union[Dict, List], obj=None, ignore=["self", "kwargs"]):
        self._obj = obj

        if (
            (not self._obj is None)
            and (not isinstance(kwargs, dict))
            and (not isinstance(kwargs, list))
        ):
            kwargs = self._obj.__dir__()

        self._ignore = ignore
        self.__parse__(kwargs, parse_ignore=["self"])

        return self.query()


# -- API-facing function: ----------------------------------------------------------------
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

    kwarg_extractor = KwargExtractor(func=func)
    return kwarg_extractor(kwargs=kwargs, obj=obj, ignore=ignore)
