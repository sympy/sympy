from sympy import cos, sin
from sympy.abc import x
from sympy.operator import DerivatedOperator

def test_DerivatedOperator():
    assert DerivatedOperator(sin)(x).doit() == cos(x)
    assert DerivatedOperator(sin)(1).doit() == cos(1)
