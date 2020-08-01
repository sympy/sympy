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
        ).rewrite(log).simplify() == 0

    one_third_e = 1*exp(1)/3
    two_thirds_e = 2*exp(1)/3
    logThirdE = log(one_third_e)
    logTwoThirdsE = log(two_thirds_e)
    lae_sum_to_e = logaddexp(logThirdE, logTwoThirdsE)
    assert lae_sum_to_e.rewrite(log) == 1
    assert lae_sum_to_e.simplify() == 1
    assert logaddexp(2, 3).simplify() == logaddexp(2, 3)  # cannot simplify with 2, 3
