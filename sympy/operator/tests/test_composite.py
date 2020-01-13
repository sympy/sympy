from sympy import Add, cos, Function, sin, symbols
from sympy.abc import x,y,z
from sympy.operator import CompositeOperator

f,g,h = symbols('f, g, h', cls=Function)

def test_CompositeOperator():
    assert CompositeOperator(sin, cos)(x) == sin(cos(x))
    assert CompositeOperator(sin, cos)(x,y) == sin(cos(x))
    assert CompositeOperator(Add, sin,cos)(x) == sin(x)+cos(x)
