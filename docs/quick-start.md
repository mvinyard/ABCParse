---
description: This guide offers a brief overview of functionality
---

# Quick Start

## Install the library

The best way to interact with our API is to use one of our official libraries:

{% tabs %}
{% tab title="pip" %}
```shell
# Install via pip
pip install adata_query
```
{% endtab %}

{% tab title="GitHub (Developer version)" %}
```sh
# Install the developer version via GitHub
git clone https://github.com/mvinyard/AnnDataQuery.git; cd ./AnnDataQuery;
pip install -e .
```
{% endtab %}
{% endtabs %}

## AnnData

This package is downstream of data loading and assumes a generally typical implementation of `adata` created using the `AnnData` package or `Scanpy`.&#x20;

```python
import anndata

h5ad_path = "/path/to/your/adata.h5ad"

adata = anndata.read_h5ad(h5ad_path)
```

Once you have some data, you are ready to interface with <mark style="color:blue;">**`adata_query`**</mark>.

## **`adata_query.`**<mark style="color:blue;">**`fetch`**</mark>:

This is probably the most useful function in the library and relies on the two functions, below. In short, this function takes a string and returns a matrix by the string, from `adata`. You can do this in grouped fashion, based on `pd.groupby`&#x20;

{% tabs %}
{% tab title="Ungrouped" %}
```python
import adata_query

key = "X_pca" # stored in adata.obsm

data = adata_query.fetch(adata = adata, key = "X_pca")
```
{% endtab %}

{% tab title="Grouped" %}
```python
import adata_query

key = "X_pca" # stored in adata.obsm
groupby = "cluster" # cell annotation in adata.obs

data = adata_query.fetch(
    adata = adata,
    key = key,
    groupby = groupby,
)
```

In this example, the returned **`data`** is now of type: <mark style="color:green;">**`List`**</mark>.
{% endtab %}
{% endtabs %}

## `adata_query.`<mark style="color:blue;">`format_data`</mark>

These functions seem trivial, but they become useful for adding flexibility into more complex workflows.&#x20;

{% tabs %}
{% tab title="numpy -> numpy" %}
For some **`data`** stored as `np.ndarray`.

```python
import adata_query

data = adata_query.format(data) # returns np.ndarray
```
{% endtab %}

{% tab title="numpy -> torch (cpu)" %}
For some **`data`** stored as `np.ndarray`.

```python
import adata_query

data = adata_query.format(data, torch = True, device = "cpu") # torch.Tensor on cpu
```
{% endtab %}

{% tab title="numpy -> torch (gpu)" %}
For some **`data`** stored as `np.ndarray`.

```python
import adata_query

data = adata_query.format(data, torch = True) # torch.Tensor on gpu, if available

# torch.Tensor can also be explicitly declared to a specific device
data = adata_query.format(data, torch = True, device = "cuda:0")

# Apple Silicon also works and will be automatically detected
data = adata_query.format(data, torch = True, device = "mps:0")
```
{% endtab %}
{% endtabs %}

## `adata_query.`<mark style="color:blue;">`locate`</mark>

I don't anticipate this function to be widely used beyond its implementation in <mark style="color:blue;">`adata_query.fetch`</mark>.

```python
import adata_query

key = "X_pca"

attr_key = adata_query.locate(adata, key = key) # attr_key = "obsm"
```

## Example notebook

Try some examples in Google Colab:

<table data-view="cards"><thead><tr><th></th><th></th><th></th><th data-hidden data-card-target data-type="content-ref"></th><th data-hidden data-card-cover data-type="files"></th></tr></thead><tbody><tr><td></td><td></td><td>Overview of common use-cases</td><td><a href="https://colab.research.google.com/github/mvinyard/AnnDataQuery/blob/main/notebooks/anndata_query_tutorial.ipynb">https://colab.research.google.com/github/mvinyard/AnnDataQuery/blob/main/notebooks/anndata_query_tutorial.ipynb</a></td><td><a href=".gitbook/assets/anndata_query_example_nb.png">anndata_query_example_nb.png</a></td></tr></tbody></table>
