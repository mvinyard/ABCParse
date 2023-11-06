
# -- set typing: --------------------------------------------------------------
from typing import Any, List, Optional, Union


# -- controller class: --------------------------------------------------------
class AsList(object):
    """Enables flexible inputs as list with type-checking."""

    def __init__(self, *args, **kwargs):
        
        """ """
        ...

    @property
    def is_list(self) -> bool:
        return isinstance(self._input, List)
    
    @property
    def _MULTIPLE_TARGET_TYPES(self):
        return isinstance(self._target_type, List)

    def _is_target_type(self, value) -> bool:
        if self._MULTIPLE_TARGET_TYPES:
            return any([isinstance(value, target_type) for target_type in self._target_type])            
        return isinstance(value, self._target_type)
    
    def _as_list(self) -> List:
        if not self.is_list:
            return [self._input]
        return self._input
    
    @property
    def list_values(self) -> List:
        return self._as_list()

    @property
    def validated_target_types(self) -> bool:
        return all([self._is_target_type(val) for val in self.list_values])

    def __call__(
        self,
        input: Union[List[Any], Any],
        target_type: Optional[Union[type, List[type]]] = None,
        *args,
        **kwargs,
    ):
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
        
        if not self._target_type is None:
            assert self.validated_target_types, "Not all values match the target type"

        return self.list_values


# -- API-facing function: -----------------------------------------------------
def as_list(
    input: Union[List[Any], Any], target_type: Optional[Union[type, List[type]]] = None, *args, **kwargs,
):
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
    _as_list = AsList()
    return _as_list(input=input, target_type=target_type)
