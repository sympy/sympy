#
# test_linsolve.py
#
#    Test the internal implementation of linsolve.
#

from sympy.testing.pytest import raises

from sympy import S, Eq, I
from sympy.abc import x, y, z

from sympy.polys.matrices.linsolve import _linsolve
from sympy.polys.solvers import PolyNonlinearError


def test__linsolve():
    assert _linsolve([], [x]) == {x:x}
    assert _linsolve([S.Zero], [x]) == {x:x}
    assert _linsolve([x-1,x-2], [x]) is None
    assert _linsolve([x-1], [x]) == {x:1}
    assert _linsolve([x-1, y], [x, y]) == {x:1, y:S.Zero}
    assert _linsolve([2*I], [x]) is None
    raises(PolyNonlinearError, lambda: _linsolve([x*(1 + x)], [x]))


def test__linsolve_float():
    eqs = [
        y - x,
        y - 0.0216 * x
    ]
    sol = {x:0, y:0}
    assert _linsolve(eqs, (x, y)) == sol

    eqs = [
        0.8*x +         0.8*z + 0.2,
        0.9*x + 0.7*y + 0.2*z + 0.9,
        0.7*x + 0.2*y + 0.2*z + 0.5
    ]
    #sol = {x:-0.69047619047619047, y:-0.52380952380952395, z:0.44047619047619047}
    sol = {x:-0.69047619047619069, y:-0.52380952380952361, z:0.44047619047619069}
    assert _linsolve(eqs, [x,y,z]) == sol

    eqs = [
        0.9*x + 0.3*y + 0.4*z + 0.6,
        0.6*x + 0.9*y + 0.1*z + 0.7,
        0.4*x + 0.6*y + 0.9*z + 0.5
    ]
    #sol = {x:-0.50285714285714289, y:-0.43809523809523804, z:-0.040000000000000063}
    sol = {x:-0.50285714285714278, y:-0.43809523809523798, z:-0.040000000000000022}
    assert _linsolve(eqs, [x,y,z]) == sol


def test__linsolve_deprecated():
    assert _linsolve([Eq(x**2, x**2+y)], [x, y]) == {x:x, y:S.Zero}
    assert _linsolve([(x+y)**2-x**2], [x]) == {x:-y/2}
    assert _linsolve([Eq((x+y)**2, x**2)], [x]) == {x:-y/2}
