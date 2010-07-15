"""Tests for algorithms for partial fraction decomposition of rational functions. """

from sympy.polys.partfrac import apart_undetermined_coeffs
from sympy.polys.polytools import Poly

from sympy.abc import x, a, b

def test_apart_undetermined_coeffs():
    p = Poly(x**2 + 1)
    q = Poly(x + 1)
    r = 2/(x + 1) + x - 1

    assert apart_undetermined_coeffs(p, q) == r

    p = Poly(2*x - 3)
    q = Poly(x**9 - x**8 - x**6 + x**5 - 2*x**2 + 3*x - 1)
    r = (-x**7 - x**6 - x**5 + 4)/(x**8 - x**5 - 2*x + 1) + 1/(x - 1)

    assert apart_undetermined_coeffs(p, q) == r

    p = Poly(1, x, domain='ZZ[a,b]')
    q = Poly((x + a)*(x + b), x, domain='ZZ[a,b]')
    r = 1/((x + b)*(a - b)) + 1/((x + a)*(b - a))

    assert apart_undetermined_coeffs(p, q) == r
