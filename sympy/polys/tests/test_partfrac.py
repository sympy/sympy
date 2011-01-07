"""Tests for algorithms for partial fraction decomposition of rational functions. """

from sympy.polys.partfrac import (
    apart_undetermined_coeffs,
    apart_full_decomposition,
    apart,
)

from sympy import S, Poly, E, pi, Matrix, Eq
from sympy.utilities.pytest import raises
from sympy.abc import x, y, a, b, c

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

    f = a*x**4 + (2*b + 2*a*c)*x**3 + (4*b*c - a**2 + a*c**2)*x**2 + (-2*a*b + 2*b*c**2)*x - b**2
    g = a**2*x**4 + (2*a*b + 2*c*a**2)*x**3 + (4*a*b*c + b**2 + a**2*c**2)*x**2 + (2*c*b**2 + 2*a*b*c**2)*x + b**2*c**2

    assert apart(f/g, x) == 1/a - 1/(x + c)**2 - b**2/(a*(a*x + b)**2)

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

    assert apart_full_decomposition(p, q) == \
        (S(1)/5)*((-x**3 + 2*x**2 - 3*x + 4)/(x**4 - x**3 + x**2 - x + 1)) + (S(1)/5)/(x + 1)
