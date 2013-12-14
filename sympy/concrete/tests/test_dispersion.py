from sympy.core import Symbol, S, oo
from sympy.concrete.dispersion import *


def test_dispersion():
    x = Symbol("x")

    fp = S(0).as_poly(x)
    assert sorted(dispersionset(fp)) == [0]

    fp = S(2).as_poly(x)
    assert sorted(dispersionset(fp)) == [0]

    fp = (x + 1).as_poly(x)
    assert sorted(dispersionset(fp)) == [0]
    assert dispersion(fp) == 0
    fp = (x*(x + 3)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 3]
    assert dispersion(fp) == 3
    fp = ((x - 3)*(x + 3)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 6]
    assert dispersion(fp) == 6
    fp = ((x + 1)*(x + 2)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 1]
    assert dispersion(fp) == 1

    fp = (x**4 - 3*x**2 + 1).as_poly(x)
    gp = fp.shift(-3)
    assert sorted(dispersionset(fp, gp)) == [2, 3, 4]
    assert dispersion(fp, gp) == 4
    assert sorted(dispersionset(gp, fp)) == []
    assert dispersion(gp, fp) == -oo

    a = Symbol("a")
    fp = (x*(3*x**2+a)*(x-2536)*(x**3+a)).as_poly(x)
    gp = fp.as_expr().subs(x, x-345).as_poly(x)
    assert sorted(dispersionset(fp, gp)) == [345, 2881]
