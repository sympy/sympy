#
# test_linsolve.py
#
#    Test the internal implementation of linsolve.
#

from sympy.testing.pytest import raises

from sympy import S, Eq, I
from sympy.abc import x, y

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


def test__linsolve_deprecated():
    assert _linsolve([Eq(x**2, x**2+y)], [x, y]) == {x:x, y:S.Zero}
    assert _linsolve([(x+y)**2-x**2], [x]) == {x:-y/2}
    assert _linsolve([Eq((x+y)**2, x**2)], [x]) == {x:-y/2}
