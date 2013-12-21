from sympy.core import Symbol, S, oo
from sympy.concrete.dispersion import *


def test_dispersion():
    x = Symbol("x")
    a = Symbol("a")

    fp = S(0).as_poly(x)
    assert sorted(dispersionset(fp)) == [0]

    fp = S(2).as_poly(x)
    assert sorted(dispersionset(fp)) == [0]

    fp = (x + 1).as_poly(x)
    assert sorted(dispersionset(fp)) == [0]
    assert dispersion(fp) == 0

    fp = ((x + 1)*(x + 2)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 1]
    assert dispersion(fp) == 1

    fp = (x*(x + 3)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 3]
    assert dispersion(fp) == 3

    fp = ((x - 3)*(x + 3)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 6]
    assert dispersion(fp) == 6

    fp = (x**4 - 3*x**2 + 1).as_poly(x)
    gp = fp.shift(-3)
    assert sorted(dispersionset(fp, gp)) == [2, 3, 4]
    assert dispersion(fp, gp) == 4
    assert sorted(dispersionset(gp, fp)) == []
    assert dispersion(gp, fp) == -oo

    fp = (x*(3*x**2+a)*(x-2536)*(x**3+a)).as_poly(x)
    gp = fp.as_expr().subs(x, x-345).as_poly(x)
    assert sorted(dispersionset(fp, gp)) == [345, 2881]
    assert sorted(dispersionset(gp, fp)) == [2191]

    gp = ((x-2)**2*(x-3)**3*(x-5)**3).as_poly(x)
    assert sorted(dispersionset(gp)) == [0, 1, 2, 3]
    assert sorted(dispersionset(gp, (gp+4)**2)) == [1, 2]

    fp = (x*(x+2)*(x-1)).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 1, 2, 3]

    # There are some difficulties if we compute over Z[a]
    # and alpha happenes to lie in Z[a] instead of simply Z.
    # Hence we can not decide if alpha is indeed integral
    # in general.

    fp = (4*x**4 + (4*a + 8)*x**3 + (a**2 + 6*a + 4)*x**2 + (a**2 + 2*a)*x).as_poly(x)
    assert sorted(dispersionset(fp)) == [0, 1]

    # For any specific value of a, the dispersion is 3*a
    # but the algorithm can not find this in general.
    # This is the point where the resultant based Ansatz
    # is superior to the current one.
    fp = (a**2*x**3 + (a**3 + a**2 + a + 1)*x).as_poly(x)
    gp = fp.as_expr().subs(x, x - 3*a).as_poly(x)
    assert sorted(dispersionset(fp, gp)) == []

    fpa = fp.as_expr().subs(a, 2).as_poly(x)
    gpa = gp.as_expr().subs(a, 2).as_poly(x)
    assert sorted(dispersionset(fpa, gpa)) == [6]
