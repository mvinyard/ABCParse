---
description: Operational class powering the locate function.
---

# AnnDataLocator

## [<mark style="color:blue;">`AnnDataLocator`</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_locator.py#L12-L138)

{% code lineNumbers="true" %}
```python
class AnnDataLocator(ABCParse.ABCParse):
    """Query available key values of AnnData. Operational class powering the `locate` function."""
    def __init__(self, *args, **kwargs) -> None:
        """
        Parameters
        ----------
        
        Returns
        -------
        None, initializes class object.
        """
        
        self._ATTRS = {}

    def _stash(self, attr: str, attr_val: np.ndarray) -> None:
        """
        Parameters
        ----------
        attr: str
        
        attr_val: np.ndarray
        
        Returns
        -------
        None, updates `self._ATTRS` and sets the (attr, attr_val) key, value pair.
        """
        self._ATTRS[attr] = attr_val
        setattr(self, attr, attr_val)

    def _intake(self, adata: anndata.AnnData) -> None:
        """
        Parameters
        ----------
        adata
        
        Returns
        -------
        
        """
        for attr in adata.__dir__():
            if "key" in attr:
                attr_val = getattr(adata, attr)()
                self._stash(attr, attr_val)
            if attr == "layers":
                attr_val = list(getattr(adata, attr))
                self._stash(attr, attr_val)

    def _cross_reference(self, passed_key: str) -> List[str]:
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        return [key for key, val in self._ATTRS.items() if passed_key in val]

    def _query_str_vals(self, query_result: List[str]) -> str:
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        return ", ".join(query_result)

    def _format_error_msg(self, key: str, query_result: List[str]) -> str:
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        if len(query_result) > 1:
            return f"Found more than one match: [{self._query_str_vals(query_result)}]"
        return f"{key} NOT FOUND"

    def _format_output_str(self, query_result: List[str]):
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        return query_result[0].split("_keys")[0]

    def _forward(self, adata: anndata.AnnData, key: str) -> str:

        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        self._intake(adata)
        query_result = self._cross_reference(passed_key=key)

        if len(query_result) != 1:
            raise KeyError(self._format_error_msg(key, query_result))

        return self._format_output_str(query_result)

    def __call__(self, adata: anndata.AnnData, key: str) -> str:
        
        """
        Parameters
        ----------
        adata: anndata.AnnData
        
        key: str
        
        Returns
        -------
        attr: str
        """

        return self._forward(adata, key)
```
{% endcode %}

**GitHub**: [<mark style="color:blue;">GitHub.com/mvinyard/AnnDataQuery/adata\_query/\_core/\_locator.py</mark>](https://github.com/mvinyard/AnnDataQuery/blob/05af47d38a677dfafa9dde41caa4ddeb479e14ba/adata\_query/\_core/\_locator.py#L13C1-L138C41)
