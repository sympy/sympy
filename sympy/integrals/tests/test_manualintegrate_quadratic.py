from sympy import symbols, sqrt, atan, log, Abs, simplify, diff
from sympy.integrals.manualintegrate import manualintegrate

def test_quadratic_denom_rule_piecewise():
    a, x = symbols('a x', real=True)
    res = manualintegrate(1 / (a + x**2), x).doit()

    # Test for a > 0 -> atan(x/sqrt(a))/sqrt(a)
    expr_pos = res.subs(a, 2)
    assert simplify(diff(expr_pos, x) - 1/(2 + x**2)) == 0

    # Test for a < 0 -> log(Abs((x - sqrt(-a))/(x + sqrt(-a))))/(2*sqrt(-a))
    expr_neg = res.subs(a, -3)
    assert simplify(diff(expr_neg, x) - 1/(-3 + x**2)) == 0

    # Test for a = 0 -> -1/x
    res_zero = manualintegrate(1 / (x**2), x).doit()
    assert simplify(res_zero + 1/x) == 0
