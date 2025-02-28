__module_name__ = "_abc_parse.py"
__doc__ = """Better abstract base class for auto-parsing inputs."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["mvinyard.ai@gmail.com"])


# -- import packages: ---------------------------------------------------------
import abc

# -- set type hints: ----------------------------------------------------------
from typing import Any, Dict, List, Optional, Tuple

# -- import internal modules: -------------------------------------------------
from . import logging

# -- Controller class: --------------------------------------------------------
class ABCParse(abc.ABC):
    _BUILT = False

    def __init__(self, *args, **kwargs) -> None:
        """
        we avoid defining things in __init__ because this subsequently
        mandates the use of `super().__init__()`
        
        Example
        -------
        ```
        class DataConfiguration(utils.ABCParse):
            def __init__(self, x=2, y=3, *args, **kwargs):
                self.__parse__(locals(), public=[None])

            def __call__(self, x=4, y=5, z=3, *args, **kwargs):
                self.__update__(locals(), private=[None])
        
        
        dc = DataConfiguration(alpha=0.2)
        dc._PARAMS
        dc(beta=0.4)
        dc._PARAMS
        dc._kwargs
        dc._PARAMS
        ```
        """
        ...
        
    def _initialize_logger(self) -> None:
        # Initialize logger with class name and logging parameters
        self._logger = logging.get_logger(
            name=self.__class__.__name__,
            level="warning",
            file_path="abcparse.log"
        )
        self._logger.debug(f"Initializing {self.__class__.__name__}")

    def __build__(self) -> None:
        self._PARAMS = {}
        self._IGNORE = ["self", "__class__"]
        self._stored_private = []
        self._stored_public = []
        self._initialize_logger()
        self._BUILT = True
        self._logger.debug("Built internal structures")

    def __set__(
        self, key: str, val: Any, public: List = [], private: List = []
    ) -> None:
        self._PARAMS[key] = val
        
        if (key in private) and (not key in public):
            self._stored_private.append(key)
            key = f"_{key}"
            self._logger.debug(f"Setting private attribute: {key}")
        else:
            self._stored_public.append(key)
            self._logger.debug(f"Setting public attribute: {key}")
        setattr(self, key, val)

    def __set_existing__(self, key: str, val: Any) -> None:
        
        passed_key = key

        if key in self._stored_private:
            key = f"_{key}"

        if passed_key == "kwargs":
            attr = getattr(self, key)
            attr.update(val)
            setattr(self, key, attr)
            self._PARAMS.update(val)
            self._logger.debug(f"Updated kwargs: {val}")
            
        elif passed_key == "args":
            attr = getattr(self, key)
            attr += val
            setattr(self, key, attr)
            self._PARAMS[passed_key] += val
            self._logger.debug(f"Updated args: {val}")
            
        else:
            self._PARAMS[passed_key] = val
            setattr(self, key, val)
            self._logger.debug(f"Updated attribute {key}: {val}")

    @property
    def _STORED(self) -> List:
        return self._stored_private + self._stored_public

    def __setup_inputs__(self, kwargs, public, private, ignore) -> Tuple[List]:
        if not self._BUILT:
            self.__build__()

        self._IGNORE += ignore
        self._logger.debug(f"Setup inputs with ignore list: {self._IGNORE}")

        if len(public) > 0:
            private = list(kwargs.keys())
            self._logger.debug(f"Public attributes specified, setting all others as private")

        return public, private

    def __parse__(
        self,
        kwargs: Dict,
        public: Optional[List] = [None],
        private: Optional[List] = [],
        ignore: Optional[List] = [],
    ) -> None:
        """
        Made to be called during `cls.__init__` of the inherited class.
        Central function of this autoparsing base class.
        
        Parameters
        ----------
        kwargs: Dict,
        public: Optional[List] = [None],
        private: Optional[List] = [],
        ignore: Optional[List] = []
        
        Returns
        -------
        None
        """
        self._logger.debug(f"Parsing kwargs: {kwargs}")
        public, private = self.__setup_inputs__(kwargs, public, private, ignore)

        for key, val in kwargs.items():
            if not key in self._IGNORE:
                self.__set__(key, val, public, private)
        
        self._logger.info(f"Parsed {len(self._PARAMS)} parameters")
        
    def __update__(
        self,
        kwargs: Dict,
        public: Optional[List] = [None],
        private: Optional[List] = [],
        ignore: Optional[List] = [],
    ) -> None:
        """
        
        To be called after __parse__ has already been called (e.g., 
        during `cls.__call__`) of the inherited class.
        
        Parameters
        ----------
        kwargs: Dict
            Typically, `locals()`
        
        public: Optional[List] = [None]
        
        private: Optional[List] = []
        
        ignore: Optional[List] = []
        
        Second-most central function of this autoparsing base class.
        
        Returns
        -------
        None
        """
        self._logger.debug(f"Updating with kwargs: {kwargs}")
        public, private = self.__setup_inputs__(kwargs, public, private, ignore)

        updated_count = 0
        new_count = 0

        for key, val in kwargs.items():
            if not (val is None) and (key in self._STORED):
                self.__set_existing__(key, val)
                updated_count += 1

            elif not (val is None) and not (key in self._IGNORE):
                self.__set__(key, val, public, private)
                new_count += 1
        
        self._logger.info(f"Updated {updated_count} existing parameters and added {new_count} new parameters")
