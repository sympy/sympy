from __future__ import annotations
from sympy.core.sympify import _sympify
from sympy.polys.polytools import Poly

from .matexpr import MatrixExpr


class CompanionMatrix(MatrixExpr):
    """A symbolic companion matrix of a polynomial.

    Examples
    ========

    >>> from sympy import Poly, Symbol, symbols
    >>> from sympy.matrices.expressions import CompanionMatrix
    >>> x = Symbol('x')
    >>> c0, c1, c2, c3, c4 = symbols('c0:5')
    >>> p = Poly(c0 + c1*x + c2*x**2 + c3*x**3 + c4*x**4 + x**5, x)
    >>> CompanionMatrix(p)
    CompanionMatrix(Poly(x**5 + c4*x**4 + c3*x**3 + c2*x**2 + c1*x + c0,
    x, domain='ZZ[c0,c1,c2,c3,c4]'))
    """
    def __new__(cls, poly):
        poly = _sympify(poly)
        if not isinstance(poly, Poly):
            raise ValueError("{} must be a Poly instance.".format(poly))
        if not poly.is_monic:
            raise ValueError("{} must be a monic polynomial.".format(poly))
        if not poly.is_univariate:
            raise ValueError(
                "{} must be a univariate polynomial.".format(poly))
        if not poly.degree() >= 1:
            raise ValueError(
                "{} must have degree not less than 1.".format(poly))

        return super().__new__(cls, poly)


    @property
    def shape(self):
        poly = self.args[0]
        size = poly.degree()
        return size, size


    def _entry(self, i, j, **kwargs):
        from sympy.core.relational import Eq
        from sympy.functions.special.tensor_functions import KroneckerDelta
        from sympy.functions.elementary.piecewise import Piecewise
        from sympy.core.numbers import Integer

        poly = self.args[0]
        n = self.cols

        if isinstance(j, (int, Integer)) and j == n - 1:
            if isinstance(i, (int, Integer)):
                return -poly.all_coeffs()[-1 - int(i)]
            return Piecewise(*[(-poly.all_coeffs()[-1 - k], Eq(i, k)) for k in range(n)])
        elif isinstance(j, (int, Integer)) and j < n - 1:
            return KroneckerDelta(i, j + 1)
        else:
            last_col = Piecewise(*[(-poly.all_coeffs()[-1 - k], Eq(i, k)) for k in range(n)])
            return Piecewise((last_col, Eq(j, n - 1)), (KroneckerDelta(i, j + 1), True))


    def as_explicit(self):
        from sympy.matrices.immutable import ImmutableDenseMatrix
        return ImmutableDenseMatrix.companion(self.args[0])
