from __future__ import print_function, division

from sympy.core.sympify import _sympify
from sympy.matrices.expressions import MatrixExpr
from sympy import S, I, sqrt, exp

class DFT(MatrixExpr):
    """ Discrete Fourier Transform

    DFT creates the matrix associated with the Discrete Fourier Transform of
    size n in the standard basis.

    Examples
    ========
    >>> from sympy.core.symbol import symbols
    >>> from sympy.matrices.expressions.fourier import DFT
    >>> M = DFT(4)
    >>> N.as_mutable().shape
    (4, 4)
    >>> M.as_mutable()
    Matrix([
    [1/2,  1/2,  1/2,  1/2],
    [1/2, -I/2, -1/2,  I/2],
    [1/2, -1/2,  1/2, -1/2],
    [1/2,  I/2, -1/2, -I/2]])
    >>> N, i, j = symbols('N, i, j')
    >>> DFT(N)[i, j] == DFT(N)[j, i]
    True
    See also:
    ============

    sympy.discrete
    https://en.wikipedia.org/wiki/DFT_matrix
"""
    def __new__(cls, n):
        n = _sympify(n)
        cls._check_dim(n)

        obj = super(DFT, cls).__new__(cls, n)
        return obj


    n = property(lambda self: self.args[0])  # type: ignore
    shape = property(lambda self: (self.n, self.n))

    def _entry(self, i, j, **kwargs):
        w = exp(-2*S.Pi*I/self.n)
        return w**(i*j) / sqrt(self.n)

    def _eval_inverse(self):
        return IDFT(self.n)

class IDFT(DFT):
    """ Inverse Discrete Fourier Transform
    IDFT creates the inverse matrix from the above.

    Examples
    ========

    >>> from sympy.matrices.expressions.fourier import DFT
    >>> Inv = IDFT(4)
    >>> Inv.as_mutable()
    Matrix([
    [1/2,  1/2,  1/2,  1/2],
    [1/2,  I/2, -1/2, -I/2],
    [1/2, -1/2,  1/2, -1/2],
    [1/2, -I/2, -1/2,  I/2]])
    """

    def _entry(self, i, j, **kwargs):
        w = exp(-2*S.Pi*I/self.n)
        return w**(-i*j) / sqrt(self.n)

    def _eval_inverse(self):
        return DFT(self.n)
