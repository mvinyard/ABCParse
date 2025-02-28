# -- set type hints: ----------------------------------------------------------
from typing import Any, List, Optional, Union

# -- import internal modules: -------------------------------------------------
from . import logging

# -- module logger: -----------------------------------------------------------
_logger = logging.get_logger(name="as_list")

# -- controller class: --------------------------------------------------------
class AsList(object):
    """Enables flexible inputs as list with type-checking."""

    def __init__(self, *args, **kwargs) -> None:
        """
        Parameters
        ----------
        *args
            type: Any

        **kwargs
        """
        self._logger = logging.get_logger(name="AsList", level="warning")
        self._logger.debug("Initializing AsList")

    @property
    def is_list(self) -> bool:
        return isinstance(self._input, List)
    
    @property
    def _MULTIPLE_TARGET_TYPES(self) -> bool:
        return isinstance(self._target_type, List)

    def _is_target_type(self, value) -> bool:
        if self._MULTIPLE_TARGET_TYPES:
            result = any([isinstance(value, target_type) for target_type in self._target_type])
            if not result:
                self._logger.debug(f"Value {value} does not match any of the target types: {self._target_type}")
            return result
        result = isinstance(value, self._target_type)
        if not result:
            self._logger.debug(f"Value {value} is not of target type: {self._target_type}")
        return result
    
    def _as_list(self) -> List[Any]:
        if not self.is_list:
            self._logger.debug(f"Converting single value to list: {self._input}")
            return [self._input]
        return self._input
    
    @property
    def list_values(self) -> List[Any]:
        return self._as_list()

    @property
    def validated_target_types(self) -> bool:
        result = all([self._is_target_type(val) for val in self.list_values])
        if result:
            self._logger.debug("All values match target type(s)")
        else:
            self._logger.warning("Not all values match target type(s)")
        return result

    def __call__(
        self,
        input: Union[List[Any], Any],
        target_type: Optional[Union[type, List[type]]] = None,
        *args,
        **kwargs,
    ) -> List[Any]:
        """
        Parameters
        ----------
        input: Union[List[Any], Any]

        target_type: Optional[Union[type, List[type]]], default = None

        Returns
        -------
        List[Any]
        """
        
        self._input = input
        self._target_type = target_type
        
        self._logger.debug(f"Processing input: {input}, target_type: {target_type}")
        
        if not self._target_type is None:
            assert self.validated_target_types, "Not all values match the target type"
            self._logger.info(f"Validated {len(self.list_values)} values against target type(s)")

        return self.list_values


# -- API-facing function: -----------------------------------------------------
def as_list(
    input: Union[List[Any], Any],
    target_type: Optional[Union[type, List[type]]] = None,
    *args,
    **kwargs,
) -> List[Any]:
    """
    Pass input to type-consistent list.
    
    Parameters
    ----------
    input: Union[List[Any], Any]

    target_type: Optional[Union[type, List[type]]], default = None
        If not all values match the target type, an AssertionError is raised.

    Returns
    -------
    List[Any]
    """
    _logger.debug(f"as_list called with input: {input}, target_type: {target_type}")
    _as_list = AsList()
    result = _as_list(input=input, target_type=target_type, *args, **kwargs)
    _logger.info(f"Converted to list with {len(result)} elements")
    return result
