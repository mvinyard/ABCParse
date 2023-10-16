---
description: Operational class powering the fetch function.
---

# AnnDataFetcher

## [<mark style="color:blue;">`AnnDataFetcher`</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_fetcher.py#L19-L60)

```python
class AnnDataFetcher(ABCParse.ABCParse):
    """Operational class powering the fetch function."""
    def __init__(self, *args, **kwargs):

        self.__parse__(locals(), public=[None])

    @property
    def _GROUPED(self):
        return self._adata.obs.groupby(self._groupby)

    def _forward(self, adata, key):
        data = getattr(adata, locate(adata, key))[key]
        return format_data(data=data, torch = self._torch, device = self._device)

    def _grouped_subroutine(self, adata, key):
        if self._as_dict:
            for group, group_df in self._GROUPED:
                yield group, self._forward(adata[group_df.index], key)
        else:
            for group, group_df in self._GROUPED:
                yield self._forward(adata[group_df.index], key)

    def __call__(
        self,
        adata: anndata.AnnData,
        key: str,
        groupby: Optional[str] = None,
        torch: bool = False,
        device: _torch.device = autodevice.AutoDevice(),
        as_dict: bool = True,
    ):
        """
        adata: anndata.AnnData [ required ]
            Annotated single-cell data object.
        
        key: str [ required ]
            Key to access a matrix in adata. For example, if you wanted to access
            adata.obsm['X_pca'], you would pass: "X_pca".
        
        groupby: Optional[str], default = None
            Optionally, one may choose to group data according to a cell-specific
            annotation in adata.obs. This would invoke returning data as List
            
        torch: bool, default = False
            Boolean indicator of whether data should be formatted as torch.Tensor. If
            False (default), data is formatted as np.ndarray.device (torch.device) =
            autodevice.AutoDevice(). Should torch=True, the device ("cpu", "cuda:N", 
            "mps:N") may be set. The default value, autodevice.AutoDevice() will 
            indicate the use of GPU, if available.

        device: torch.device, default = autodevice.AutoDevice()
            
    
        as_dict: bool, default = True
            Only relevant when `groupby` is not None. Boolean indicator to return
            data in a Dict where the key for each value corresponds to the respective
            `groupby` value. If False, returns List.
        """

        self.__update__(locals(), public=[None])

        if hasattr(self, "_groupby"):
            if self._as_dict:
                return dict(self._grouped_subroutine(adata, key))
            return list(self._grouped_subroutine(adata, key))
        return self._forward(adata, key)
```

**GitHub**: [<mark style="color:blue;">GitHub.com/mvinyard/AnnDataQuery/adata\_query/\_core/\_fetcher.py</mark>](https://github.com/mvinyard/AnnDataQuery/blob/05af47d38a677dfafa9dde41caa4ddeb479e14ba/adata\_query/\_core/\_fetcher.py#L20C1-L85C41)
