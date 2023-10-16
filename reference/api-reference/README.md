---
description: Description of functions and the classes on which they are built
---

# API Reference

## Functions

Functional wrappers of the classes below. These are the main tools of the API.

<details>

<summary><a href="fetch/"><mark style="color:blue;"><code>adata_query.fetch</code></mark></a></summary>

```python
data = fetch(
    adata: anndata.AnnData,
    key: str,
    groupby: Optional[str] = None,
    torch: bool = False,
    device: torch.device = autodevice.AutoDevice(),
    *args,
    **kwargs,
)
```

</details>

<details>

<summary><a href="format_data/"><mark style="color:blue;"><code>adata_query.format_data</code></mark></a></summary>

```python
data = format_data(
    data: Union[np.ndarray, torch.Tensor],
    torch: bool = False,
    device: torch.device = autodevice.AutoDevice(),
    *args,
    **kwargs,
)
```

</details>

<details>

<summary><a href="locate/"><mark style="color:blue;"><code>adata_query.locate</code></mark></a></summary>

```python
attr_key = adata_query.locate(adata: anndata.AnnData, key: str)
```

</details>

{% hint style="info" %}
Expand each function for a quick view. Click on the function titles for full documentation.
{% endhint %}

## Classes

Everything related to the operational classes:

{% content-ref url="fetch/anndatafetcher.md" %}
[anndatafetcher.md](fetch/anndatafetcher.md)
{% endcontent-ref %}

{% content-ref url="format_data/dataformatter.md" %}
[dataformatter.md](format\_data/dataformatter.md)
{% endcontent-ref %}

{% content-ref url="locate/anndatalocator.md" %}
[anndatalocator.md](locate/anndatalocator.md)
{% endcontent-ref %}
