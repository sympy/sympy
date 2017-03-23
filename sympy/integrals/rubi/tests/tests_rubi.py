from sympy.integrals.rubi.Algebra_Integration_rules import int_alg_1_1
from sympy import S, Symbol

x = Symbol('x')
def test_int_alg_1_1():
    assert int_alg_1_1((5+6*x)**3,x) == (6*x + 5)**4/24