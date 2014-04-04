from __future__ import print_function, division

from .matexpr import MatrixExpr, ShapeError, Identity
from sympy import Pow, S, Basic
from sympy.core.sympify import _sympify


class MatPow(MatrixExpr):

    def __new__(cls, base, exp):
        base = _sympify(base)
        if not base.is_Matrix:
            raise TypeError("Function parameter should be a matrix")
        exp = _sympify(exp)
        return super(MatPow, cls).__new__(cls, base, exp)

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    @property
    def shape(self):
        return self.base.shape

    def _entry(self, i, j):
        if self.exp.is_Integer:
            # Make an explicity MatMul out of the MatPow
            return MatMul(*[self.base for k in range(self.exp)])._entry(i, j)

from .matmul import MatMul
