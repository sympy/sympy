"""Tests for the ``sympy.simplify.cse_diff.py`` module."""

import pytest

from sympy.core.symbol import Symbol
from sympy.core.numbers import Integer
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.trigonometric import sin, tan
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.physics.mechanics import dynamicsymbols
from sympy.simplify.cse_diff import _forward_jacobian
from sympy.simplify.simplify import simplify

w = Symbol('w')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

q1, q2, q3 = dynamicsymbols('q1 q2 q3')

zero = Integer(0)
one = Integer(1)
two = Integer(2)
neg_one = Integer(-1)

@pytest.mark.parametrize(
    'expr, wrt',
    [
        ([zero], [x]),
        ([one], [x]),
        ([two], [x]),
        ([neg_one], [x]),
        ([x], [x]),
        ([y], [x]),
        ([x + y], [x]),
        ([x*y], [x]),
        ([x**2], [x]),
        ([x**y], [x]),
        ([exp(x)], [x]),
        ([sin(x)], [x]),
        ([tan(x)], [x]),
        ([zero, one, x, y, x*y, x + y], [x, y]),
        ([((x/y) + sin(x/y) - exp(y))*((x/y) - exp(y))], [x, y]),
        ([w*tan(y*z)/(x - tan(y*z)), w*x*tan(y*z)/(x - tan(y*z))], [w, x, y, z]),
        ([q1**2 + q2, q2**2 + q3, q3**2 + q1], [q1, q2, q3])
    ]
)
def test_forward_jacobian(expr, wrt):
    expr = ImmutableDenseMatrix([expr]).T
    wrt = ImmutableDenseMatrix([wrt]).T
    jacobian = _forward_jacobian(expr, wrt, as_cse_expr=True)
    zeros = ImmutableDenseMatrix.zeros(*jacobian.shape)
    assert simplify(jacobian - expr.jacobian(wrt)) == zeros
