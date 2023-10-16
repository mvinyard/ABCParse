---
description: Locate the attribute data container in AnnData.
---

# locate

## [<mark style="color:blue;">`locate`</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_locator.py#L141-L154)

Given, <mark style="color:blue;">`adata`</mark> and a key that points to a specific matrix stored in <mark style="color:blue;">`adata`</mark>,  return the data, formatted either as `np.ndarray` or `torch.Tensor`. If formatted as `torch.Tensor`, `device` may be specified based on available devices.

Uses the class: [<mark style="color:blue;">**`adata_query._core.AnnDataLocator`**</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_locator.py#L13C1-L138).

{% code overflow="wrap" %}
```python
attr_key = locate(adata: anndata.AnnData, key: str)
```
{% endcode %}

<details>

<summary>Parameters</summary>

**`adata`** (<mark style="color:blue;">`anndata.AnnData`</mark>)

* Annotated single-cell data object. See [`anndata.AnnData`](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html#anndata.AnnData).

**`key`** (<mark style="color:blue;">`str`</mark>)

* Key to access a matrix in `adata`. For example, if you wanted to access `adata.obsm['X_pca']`, you would pass: `"X_pca"`.

</details>

<details>

<summary>Returns</summary>

**`attr_key`**(<mark style="color:blue;">`str`</mark>)

* Attribute of `adata` containing the passed `key`

</details>

<details>

<summary>Exceptions</summary>

There are two possible `KeyErrors` that may arise when using this function. The first is the requirement that the passed `key` is found in any attribute of `adata`. Should it not be found, the following `KeyError` will be returned:

<mark style="color:red;">**`KeyError: {key} NOT FOUND`**</mark>

The second possible `KeyError` that may arise indicates the discovery of multiple `adata` attributes containing a key-value pair where the passed key is observed in both attributes. For example, if the passed key matches a column name in `obs` and a key in `adata.obsm`, the following `KeyError` would arise:

<mark style="color:red;">**`KeyError: Found more than one match: ['obs', 'obsm']`**</mark>

In the event of the second `KeyError` in which multiple locations are found, rename the `key`-`value` pairs in `adata`, such that each is unique.

* Attribute of `adata` containing the passed `key`

</details>

**Source code**: [<mark style="color:blue;">GitHub.com/mvinyard/AnnDataQuery/adata\_query/\_core/\_locator.py</mark>](https://github.com/mvinyard/AnnDataQuery/blob/05af47d38a677dfafa9dde41caa4ddeb479e14ba/adata\_query/\_core/\_locator.py#L141-L161)
