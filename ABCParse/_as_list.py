

# -- import local dependenecies: ----------------------------------------------
from ._abc_parse import ABCParse


# -- set typing: --------------------------------------------------------------
from typing import Union, List, Any


# -- controller class: --------------------------------------------------------
class AsList(ABCParse):
    """Enables flexible inputs as list with type-checking."""

    def __init__(self, *args, **kwargs):
        self.__parse__(locals(), public=[None])

    @property
    def is_list(self) -> bool:
        return isinstance(self._input, List)

    def _is_target_type(self, value) -> bool:
        return isinstance(value, self._target_type)

    @property
    def list_values(self):
        return self._as_list()

    @property
    def validated_target_types(self):
        return all([self._is_target_type(val) for val in self.list_values])

    def _as_list(self):
        if not self.is_list:
            return [self._input]
        return self._input

    def __call__(
        self,
        input: Union[List[Any], Any],
        target_type: type = str,
        *args,
        **kwargs,
    ):
        """
        Parameters
        ----------
        input: Union[List[Any], Any]

        target_type: typel, default = str

        Returns
        -------
        List[Any]
        """
        self.__update__(locals(), public=[None])

        assert self.validated_target_types, "Not all values match the target type"

        return self.list_values


# -- API-facing function: -----------------------------------------------------
def as_list(
    input: Union[List[Any], Any], target_type: type = str, *args, **kwargs,
):
    """
    Pass input to type-consistent list.
    
    Parameters
    ----------
    input: Union[List[Any], Any]

    target_type: typel, default = str
        If not all values match the target type, an AssertionError is raised.

    Returns
    -------
    List[Any]
    """
    _as_list = AsList()
    return _as_list(input=input, target_type=target_type)
