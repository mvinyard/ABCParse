---
description: Format data to interface with numpy or torch, on a specified device.
---

# format\_data

## [<mark style="color:blue;">`format_data`</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_formatter.py#L81C1-L86C32)

Given, <mark style="color:blue;">`adata`</mark> and a key that points to a specific matrix stored in <mark style="color:blue;">`adata`</mark>,  return the data, formatted either as `np.ndarray` or `torch.Tensor`. If formatted as `torch.Tensor`, `device` may be specified based on available devices.

Uses the class: [<mark style="color:blue;">**adata\_query.\_core.DataFormatter**</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_formatter.py#L15-L77).

{% code overflow="wrap" %}
```python
data = format_data(
    data: Union[np.ndarray, torch.Tensor],
    torch: bool = False,
    device: torch.device = autodevice.AutoDevice(),
    *args,
    **kwargs,
)
```
{% endcode %}

<details>

<summary>Parameters</summary>

**`data`** (<mark style="color:blue;">`Union[np.ndarray, torch.Tensor]`</mark>)

* Input data that should be formatted. Typically an `np.ndarray`, `torch.Tensor`, or `ArrayView`.

**`torch`** (<mark style="color:blue;">`bool`</mark>) `=`` `<mark style="color:green;">**`False`**</mark>

* Boolean indicator of whether data should be formatted as `torch.Tensor`. If <mark style="color:green;">**`False`**</mark> (default), data is formatted as np.ndarray.

**`device`** (<mark style="color:blue;">`torch.device`</mark>) `= autodevice.AutoDevice()`

* Should **`torch`**`=`<mark style="color:green;">**`True`**</mark>, the device (`"cpu"`, `"cuda:N"`, `"mps:N"`) may be set. The default value, `autodevice.AutoDevice()` will indicate the use of GPU, if available.

</details>

<details>

<summary>Returns</summary>

**`data`**(<mark style="color:blue;">`Union[torch.Tensor, np.ndarray])`</mark>

* Formatted data as `np.ndarray` or `torch.Tensor`.
* If <mark style="color:blue;">**`torch`**</mark>**`=`**<mark style="color:green;">**`True`**</mark>, the <mark style="color:blue;">**`torch.Tensor`**</mark> is allocated to the device indicated by the <mark style="color:blue;">**`device`**</mark> argument.

</details>

**Source code**: [GitHub.com/mvinyard/AnnDataQuery/adata\_query/\_core/\_formatter.py](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_formatter.py#L81-L86)

