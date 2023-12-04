# Issues

## 1. 2023-10-27-23:46

```python
def calSumPsiJ(self, x, y) -> float:
    nodal_part = self.calSumPsiJNodalPart(x, y)
    linear_part = self.calSumPsiJLinearPart(x, y)
    return nodal_part + linear_part
    pass
```

linear_part 会出现nan 

fixed. 