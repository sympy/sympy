from sympy import Expr, ImmutableDenseNDimArray
from sympy.core.sympify import _sympify


class ZeroArray(Expr):
    """
    Symbolic array of zeros. Equivalent to ``ZeroMatrix`` for matrices.
    """

    def __new__(cls, *shape):
        shape = map(_sympify, shape)
        obj = Expr.__new__(cls, *shape)
        return obj

    @property
    def shape(self):
        return self._args

    def as_explicit(self):
        if any(not i.is_Integer for i in self.shape):
            raise ValueError("Cannot return explicit form for symbolic shape.")
        return ImmutableDenseNDimArray.zeros(*self.shape)
