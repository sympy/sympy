from sympy import S, symbols
from sympy.stats import E
from sympy.stats.stochastic_process_types import BernoulliProcess

def test_BernoulliProcess():
    B = BernoulliProcess('x', S(1)/3)
    x, y = symbols('x y', integer=True)
    assert E(B[1]) == S(1)/3
    assert E(B[x + y], B[x]) == S(1)/3
