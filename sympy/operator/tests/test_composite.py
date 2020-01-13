from sympy import Add, cos, sin
from sympy.abc import x,y
from sympy.operator import CompositeOperator

def test_CompositeOperator():
    assert CompositeOperator(sin, cos)(x) == sin(cos(x))
    assert CompositeOperator(sin, cos)(x,y) == sin(cos(x))
    assert CompositeOperator(Add, sin,cos)(x) == sin(x)+cos(x)
