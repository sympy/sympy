"""Tests for algorithms for partial fraction decomposition of rational functions. """

from sympy.polys.partfrac import (
    apart_undetermined_coeffs,
    apart_full_decomposition,
    apart,
)

from sympy import S, Poly, raises, E, pi, Matrix, Eq
from sympy.abc import x, y, a, b

def test_apart():
    assert apart(1) == 1
    assert apart(1, x) == 1

    f, g = (x**2 + 1)/(x + 1), 2/(x + 1) + x - 1

    assert apart(f, full=False) == g
    assert apart(f, full=True) == g

    f, g = 1/(x+2)/(x+1), 1/(1 + x) - 1/(2 + x)

    assert apart(f, full=False) == g
    assert apart(f, full=True) == g

    f, g = 1/(x+1)/(x+5), -1/(5 + x)/4 + 1/(1 + x)/4

    assert apart(f, full=False) == g
    assert apart(f, full=True) == g

    raises(NotImplementedError, "apart(1/(x + 1)/(y + 2))")

    assert apart((E*x+2)/(x-pi)*(x-1), x) in [
        2 - E + E*pi + E*x - 1/(x - pi)*( 2 - 2*pi + E*pi - E*pi**2),
        2 - E + E*pi + E*x + 1/(x - pi)*(-2 + 2*pi - E*pi + E*pi**2),
    ]

    M = Matrix(2, 2, lambda i, j: 1/(x - (i+1))/(x - (1-j)))

    assert apart(M) in [
        Matrix([
            [(x - 1)**(-2),         -1/x - 1/(1 - x)          ],
            [1/(1 - x) - 1/(2 - x), -S.Half/x - S.Half/(2 - x)],
        ]),
        Matrix([
            [(x - 1)**(-2),          -1/x + 1/(x - 1)          ],
            [-1/(x - 1) + 1/(x - 2), -S.Half/x + S.Half/(x - 2)],
        ]),
    ]

    assert apart(Eq((x**2 + 1)/(x + 1), x), x) == Eq(x - 1 + 2/(x + 1), x)

def test_apart_undetermined_coeffs():
    p = Poly(2*x - 3)
    q = Poly(x**9 - x**8 - x**6 + x**5 - 2*x**2 + 3*x - 1)
    r = (-x**7 - x**6 - x**5 + 4)/(x**8 - x**5 - 2*x + 1) + 1/(x - 1)

    assert apart_undetermined_coeffs(p, q) == r

    p = Poly(1, x, domain='ZZ[a,b]')
    q = Poly((x + a)*(x + b), x, domain='ZZ[a,b]')
    r = 1/((x + b)*(a - b)) + 1/((x + a)*(b - a))

    assert apart_undetermined_coeffs(p, q) == r

def test_apart_full_decomposition():
    p = Poly(1, x)
    q = Poly(x**5 + 1, x)

    assert str(apart_full_decomposition(p, q)) in [
        "RootSum(x**4 - x**3 + x**2 - x + 1, Lambda(_a, -1/5/(x - _a)*_a)) + 1/(5*(1 + x))",
        "RootSum(x**4 - x**3 + x**2 - x + 1, Lambda(_a, -_a/(5*(x - _a)))) + 1/(5*(1 + x))",
    ]

