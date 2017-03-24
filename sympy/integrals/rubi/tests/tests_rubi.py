from sympy.integrals.rubi.Algebra_Integration_rules import int_alg_1_1
from sympy.integrals.rubi.Algebra_Integration_rules import rubi_integrate
from sympy import S, Symbol

x = Symbol('x')
def test_int_alg_1_1():
    assert int_alg_1_1((5 + 6*x)**3, x) == (6*x + 5)**4/24
    assert rubi_integrate((5 + 6*x)**3, x) == (6*x + 5)**4/24
    assert rubi_integrate((6*x)**3, x) == 54*x**4
    assert rubi_integrate(2*(2 + 3*x)**5, x) == (3*x + 2)**6/9
