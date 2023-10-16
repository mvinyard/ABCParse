---
description: Operational class powering the format_data function.
---

# DataFormatter

## [<mark style="color:blue;">`DataFormatter`</mark>](https://github.com/mvinyard/AnnDataQuery/blob/fa2b5728c0c24752b03a1335866b576179748cbb/adata\_query/\_core/\_formatter.py#L15-L77)

{% code overflow="wrap" lineNumbers="true" %}
```python
class DataFormatter(ABCParse.ABCParse):
    """Format data to interface with numpy or torch, on a specified device."""
    def __init__(self, data: Union[_torch.Tensor, np.ndarray], *args, **kwargs):
        self.__parse__(locals())

    @property
    def device_type(self) -> str:
        """Returns device type"""
        if hasattr(self.data, "device"):
            return self.data.device.type
        return "cpu"

    @property
    def is_ArrayView(self) -> bool:
        """Checks if device is of type ArrayView"""
        return isinstance(self.data, anndata._core.views.ArrayView)

    @property
    def is_numpy_array(self) -> bool:
        """Checks if device is of type np.ndarray"""
        return isinstance(self.data, np.ndarray)

    @property
    def is_torch_Tensor(self) -> bool:
        """Checks if device is of type torch.Tensor"""
        return isinstance(self.data, _torch.Tensor)

    @property
    def on_cpu(self) -> bool:
        """Checks if device is on cuda or mps"""
        return self.device_type == "cpu"

    @property
    def on_gpu(self) -> bool:
        """Checks if device is on cuda or mps"""
        return self.device_type in ["cuda", "mps"]

    def to_numpy(self) -> np.ndarray:
        """Sends data to np.ndarray"""
        if self.is_torch_Tensor:
            if self.on_gpu:
                return self.data.detach().cpu().numpy()
            return self.data.numpy()
        elif self.is_ArrayView:
            return self.data.toarray()
        return self.data

    def to_torch(self, device=autodevice.AutoDevice()) -> _torch.Tensor:
        """
        Parameters
        ----------
        device: torch.device

        Returns
        -------
        torch.Tensor
        """
        self.__update__(locals())

        if self.is_torch_Tensor:
            return self.data.to(device)
        elif self.is_ArrayView:
            self.data = self.data.toarray()
        return _torch.Tensor(self.data).to(device)
```
{% endcode %}

**GitHub**: [<mark style="color:blue;">GitHub.com/mvinyard/AnnDataQuery/adata\_query/\_core/\_formatter.py</mark>](https://github.com/mvinyard/AnnDataQuery/blob/05af47d38a677dfafa9dde41caa4ddeb479e14ba/adata\_query/\_core/\_formatter.py#L82-L104)
