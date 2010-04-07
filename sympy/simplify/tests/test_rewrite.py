from sympy import (sin, cos, exp, cot, sqrt, S, I, E, pi, symbols, Function,
    Matrix, Eq, RootSum, Lambda)
from sympy.simplify import apart, together
from sympy.integrals import integrate
from sympy.utilities.pytest import XFAIL, raises

x,y,z,n = symbols('xyzn')

def test_has():
    assert cot(x).has(x)
    assert cot(x).has(cot)
    assert not cot(x).has(sin)
    assert sin(x).has(x)
    assert sin(x).has(sin)
    assert not sin(x).has(cot)

def test_sin_exp_rewrite():
    assert sin(x).rewrite(sin, exp) == -I/2*(exp(I*x)-exp(-I*x))
    assert sin(x).rewrite(sin, exp).rewrite(exp, sin) == sin(x)
    assert cos(x).rewrite(cos, exp).rewrite(exp, cos) == cos(x)
    assert (sin(5*y) - sin(2*x)).rewrite(sin, exp).rewrite(exp, sin) == sin(5*y) - sin(2*x)
    assert sin(x+y).rewrite(sin, exp).rewrite(exp, sin) == sin(x+y)
    assert cos(x+y).rewrite(cos, exp).rewrite(exp, cos) == cos(x+y)
    # This next test currently passes... not clear whether it should or not?
    assert cos(x).rewrite(cos, exp).rewrite(exp, sin) == cos(x)

def test_apart():
    raises(ValueError, "apart(1/(x+1)/(y+2))")

    assert apart(1) == 1
    assert apart(1, x) == 1

    assert apart(1/(x+2)/(x+1)) == 1/(1 + x) - 1/(2 + x)
    assert apart(1/(x+1)/(x+5)) == -1/(5 + x)/4 + 1/(1 + x)/4

    f = apart(1/(x-y)/(x-z), x)

    assert f.subs({y:1,z:2}) == apart(1/(x-1)/(x-2), x)

    assert apart((E*x+2)/(x-pi)*(x-1), x) in [
        2 - E + E*pi + E*x - 1/(x - pi)*( 2 - 2*pi + E*pi - E*pi**2),
        2 - E + E*pi + E*x + 1/(x - pi)*(-2 + 2*pi - E*pi + E*pi**2),
    ]

    M = Matrix(2, 2, lambda i, j: 1/(x-(i+1))/(x-(1-j)))

    assert apart(M, x) in [
            Matrix([
                [(x-1)**(-2),     -1/x-1/(1-x)          ],
                [1/(1-x)-1/(2-x), -S.Half/x-S.Half/(2-x)],
            ]),
            Matrix([
                [(-1+x)**(-2),     -1/x+1/(-1+x)          ],
                [-1/(-1+x)+1/(-2+x), -S.Half/x+S.Half/(-2+x)],
            ]),
            ]

    assert apart(Eq((x**2+1)/(x+1), sin(x)), x) == \
        Eq(x - 1 + 2/(x+1), sin(x))

    assert str(apart(1/(1+x**5), x)) in [
        "RootSum(x**4 - x**3 + x**2 - x + 1, Lambda(_a, -1/5/(x - _a)*_a)) + 1/(5*(1 + x))",
        "RootSum(x**4 - x**3 + x**2 - x + 1, Lambda(_a, -_a/(5*(x - _a)))) + 1/(5*(1 + x))",
    ]

