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

    # This should give the exact answer:
    eqs = [
        y - x,
        y - 0.0216 * x
    ]
    sol = {x:0.0, y:0.0}
    assert _linsolve(eqs, (x, y)) == sol

    # Other cases should be close to eps

    def all_close(sol1, sol2, eps=1e-15):
        close = lambda a, b: abs(a - b) < eps
        assert sol1.keys() == sol2.keys()
        return all(close(sol1[s], sol2[s]) for s in sol1)

    eqs = [
        0.8*x +         0.8*z + 0.2,
        0.9*x + 0.7*y + 0.2*z + 0.9,
        0.7*x + 0.2*y + 0.2*z + 0.5
    ]
    sol_exact = {x:-29/42, y:-11/21, z:37/84}
    sol_linsolve = _linsolve(eqs, [x,y,z])
    assert all_close(sol_exact, sol_linsolve)

    eqs = [
        0.9*x + 0.3*y + 0.4*z + 0.6,
        0.6*x + 0.9*y + 0.1*z + 0.7,
        0.4*x + 0.6*y + 0.9*z + 0.5
    ]
    sol_exact = {x:-88/175, y:-46/105, z:-1/25}
    sol_linsolve = _linsolve(eqs, [x,y,z])
    assert all_close(sol_exact, sol_linsolve)

    eqs = [
        x*(0.7 + 0.6*I) + y*(0.4 + 0.7*I) + z*(0.9 + 0.1*I) + 0.5,
        0.2*I*x + 0.2*I*y + z*(0.9 + 0.2*I) + 0.1,
        x*(0.9 + 0.7*I) + y*(0.9 + 0.7*I) + z*(0.9 + 0.4*I) + 0.4,
    ]
    sol_exact = {
        x:-6157/7995 - 411/5330*I,
        y:8519/15990 + 1784/7995*I,
        z:-34/533 + 107/1599*I,
    }
    sol_linsolve = _linsolve(eqs, [x,y,z])
    assert all_close(sol_exact, sol_linsolve)


def test__linsolve_deprecated():
    assert _linsolve([Eq(x**2, x**2+y)], [x, y]) == {x:x, y:S.Zero}
    assert _linsolve([(x+y)**2-x**2], [x]) == {x:-y/2}
    assert _linsolve([Eq((x+y)**2, x**2)], [x]) == {x:-y/2}
