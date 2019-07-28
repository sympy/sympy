from __future__ import print_function, division

from sympy.core.sympify import _sympify
from sympy.core.decorators import deprecated

from .matexpr import MatrixExpr
from .inverse import Inverse

from sympy import S, I, sqrt, exp

class DFTMatrix(MatrixExpr):
    """Discrete Fourier transformation matrix"""
    n = property(lambda self: self.args[0])
    shape = property(lambda self: (self.n, self.n))

    def __new__(cls, n, fourier_params=(0, -1)):
        n = _sympify(n)
        if n.is_number and not (n.is_integer and n.is_positive):
            raise ValueError(
                'Matrix size {} should be specified as a positive integer.')

        (a, b) = _sympify(fourier_params)

        return MatrixExpr().__new__(cls, n, a, b)

    def _entry(self, i, j):
        n = self.rows
        a, b = self.args[1], self.args[2]
        w = exp(2*S.Pi*I*b / n)

        return w**(i*j) / sqrt(n**(1-a))

    def _eval_inverse(self):
        n = self.rows
        a, b = self.args[1], self.args[2]

        if b.is_number and b.is_integer:
            if b.gcd(n) != 1:
                raise ValueError(
                    'The DFT Matrix with the specification is not invertible')
            if not b.is_integer:
                return Inverse(self)

            return DFTMatrix(n, fourier_params=(-a, -b))

@deprecated(
    issue=99999,
    useinstead="DFTMatrix",
    deprecated_since_version="1.5")
def DFT(n):
    """Discrete Fourier transformation matrix."""
    return DFTMatrix(n)


@deprecated(
    issue=99999,
    useinstead="DFTMatrix.inverse",
    deprecated_since_version="1.5")
def IDFT(n):
    return DFTMatrix(n).inverse()
