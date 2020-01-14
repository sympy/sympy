from sympy import cos, sin
from sympy.abc import x
from sympy.operator.derivative import DerivatedOp

def test_DerivatedOp():
    assert DerivatedOp(sin)(x).doit() == cos(x)
    assert DerivatedOp(sin)(1).doit() == cos(1)
