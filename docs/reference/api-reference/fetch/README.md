---
description: Fetch data from AnnData, flexibly.
---

# fetch

## [<mark style="color:blue;">`fetch`</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_fetcher.py#L62-L82C6)

Given, <mark style="color:blue;">`adata`</mark> and a key that points to a specific matrix stored in <mark style="color:blue;">`adata`</mark>,  return the data, formatted either as `np.ndarray` or `torch.Tensor`. If formatted as `torch.Tensor`, `device` may be specified based on available devices.

Uses the class: [<mark style="color:blue;">**`adata_query._core.AnnDataFetcher`**</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_fetcher.py#L19-L60https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_fetcher.py#L19-L60).

{% code overflow="wrap" %}
```python
data = fetch(
    adata: anndata.AnnData,
    key: str,
    groupby: Optional[str] = None,
    torch: bool = False,
    device: torch.device = autodevice.AutoDevice(),
    as_dict: bool = True,
    *args,
    **kwargs,
)
```
{% endcode %}

<details>

<summary>Parameters</summary>

**`adata`** (<mark style="color:blue;">`anndata.AnnData`</mark>)

* Annotated single-cell data object. See [`anndata.AnnData`](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html#anndata.AnnData)

**`key`** (<mark style="color:blue;">`str`</mark>)

* Key to access a matrix in `adata`. For example, if you wanted to access `adata.obsm['X_pca']`, you would pass: `"X_pca"`.

**`groupby`** (<mark style="color:blue;">`Optional[str]`</mark>) `=`` `<mark style="color:green;">**`None`**</mark>

* Optionally, one may choose to group data according to a cell-specific annotation in `adata.obs`. This would invoke returning data as **`List`**

**`torch`** (<mark style="color:blue;">`bool`</mark>) `=`` `<mark style="color:green;">**`False`**</mark>

* Boolean indicator of whether data should be formatted as `torch.Tensor`. If <mark style="color:green;">**`False`**</mark> (default), data is formatted as np.ndarray.

**`device`** (<mark style="color:blue;">`torch.device`</mark>) `= autodevice.AutoDevice()`

* Should **`torch`**`=`<mark style="color:green;">**`True`**</mark>, the device (`"cpu"`, `"cuda:N"`, `"mps:N"`) may be set. The default value, `autodevice.AutoDevice()` will indicate the use of GPU, if available.

**`as_dict`** (<mark style="color:blue;">`bool`</mark>) `=`<mark style="color:green;">**`True`**</mark>

* Only relevant when `groupby` is not <mark style="color:green;">**None**</mark>. Boolean indicator to return data in a **`Dict`** where the key for each value corresponds to the respective `groupby` value. If False, returns **`List`**.

</details>

<details>

<summary>Returns</summary>

**`data`**(<mark style="color:blue;">`Union[torch.Tensor, np.ndarray, List[Union[torch.Tensor, np.ndarray]], Dict[Union[str, int], Union[torch.Tensor, np.ndarray]]`</mark>)

* Formatted data as `np.ndarray` or torch.Tensor. If `torch=`<mark style="color:green;">**`True`**</mark> the `torch.Tensor` is allocated to the device indicated by the device argument. If `groupby` is passed, returned as **`Dict`**`[`**`Union`**`[str, int], np.ndarray]` or **`Dict`**`[`**`Union`**`[str, int]. torch.Tensor]`. If groupby is passed and `as_dict =`` `<mark style="color:green;">**`False`**</mark>, returns **`List`**`[np.ndarray]` or **`List`**`[tobrch.Tensor]`.

</details>

**Source code**: [GitHub.com/mvinyard/AnnDataQuery/adata\_query/\_core/\_fetcher.py](https://github.com/mvinyard/AnnDataQuery/blob/05af47d38a677dfafa9dde41caa4ddeb479e14ba/adata\_query/\_core/\_fetcher.py#L87-L152)
