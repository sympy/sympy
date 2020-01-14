from sympy import cos, sin
from sympy.abc import x,y
from sympy.operator.composite import CompositeOp

def test_CompositeOp():
    assert CompositeOp(sin, cos)(x) == sin(cos(x))
    assert CompositeOp(sin, cos)(x,y) == sin(cos(x))
