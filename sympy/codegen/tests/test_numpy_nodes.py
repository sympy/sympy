from itertools import product
from sympy import symbols, exp, log
from sympy.codegen.numpy_nodes import logaddexp

x, y, z = symbols('x y z')

def test_logaddexp():
    lae_xy = logaddexp(x, y)
    ref_xy = log(exp(x) + exp(y))
    for wrt, deriv_order in product([x, y, z], range(0, 3)):
        assert (
            lae_xy.diff(wrt, deriv_order) -
            ref_xy.diff(wrt, deriv_order)
        ).simplify() == 0
