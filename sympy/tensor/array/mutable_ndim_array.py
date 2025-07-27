from sympy.tensor.array.ndim_array import NDimArray
from typing_extensions import Self


class MutableNDimArray(NDimArray):

    def as_immutable(self):
        raise NotImplementedError("abstract method")

    def as_mutable(self) -> Self:
        return self

    def _sympy_(self):
        return self.as_immutable()
