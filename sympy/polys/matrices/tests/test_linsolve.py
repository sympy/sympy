#
# test_linsolve.py
#
#    Test the internal implementation of linsolve.
#


from sympy import S
from sympy.abc import x

from sympy.polys.matrices.linsolve import _linsolve


def test__linsolve():
    assert _linsolve([], [x]) == {x:x}
    assert _linsolve([S.Zero], [x]) == {x:x}
