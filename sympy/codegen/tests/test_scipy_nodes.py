from itertools import product
from sympy.core.symbol import symbols
from sympy.functions.elementary.trigonometric import cos
from sympy.core.numbers import pi
from sympy.codegen.scipy_nodes import cosm1

x, y, z = symbols('x y z')


def test_cosm1():
    cm1_xy = cosm1(x*y)
    ref_xy = cos(x*y) - 1
    for wrt, deriv_order in product([x, y, z], range(0, 3)):
        assert (
            cm1_xy.diff(wrt, deriv_order) -
            ref_xy.diff(wrt, deriv_order)
        ).rewrite(cos).simplify() == 0

    expr_minus2 = cosm1(pi)
    assert expr_minus2.rewrite(cos) == -2
    assert cosm1(3.14).simplify() == cosm1(3.14)  # cannot simplify with 3.14
