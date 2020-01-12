from sympy import Add, cos, Function, sin, symbols
from sympy.abc import x,y,z
from sympy.operator import Operator

f,g,h = symbols('f, g, h', cls=Function)

def test_Operator():
    assert Operator(1)(x,y,z) == 1
    assert Operator(sin)(x) == sin(x)
    assert Operator(cos)(y) == cos(y)
    assert Operator(Add)(x,y,z) == x+y+z
    assert Operator(f)(x) == f(x)
    assert Operator(g)(x,y) == g(x,y)
    assert Operator(h(x))(x,y) == h(x)

def test_AppliedOperator1():
    """
    Cases where AppliedOperator is force-evaluated.
    """
    assert Operator(1)(x,y, evaluate=False) == 1
    assert Operator(f(x))(y,z, evaluate=False) == f(x)
    assert Operator(sin)(x, evaluate=False) == sin(x)
    assert Operator(Add)(x,x, evaluate=False) == Add(x,x, evaluate=False)

def test_AppliedOperator2():
    """
    Cases where AppliedOperator is not force-evaluated.
    """
    op1 = sin+cos
    op2 = sin/cos
    assert op1(x, evaluate=False) != sin(x)+cos(x)
    assert op1(x, evaluate=False).doit() == sin(x)+cos(x)
    assert op2(x, evaluate=False) != sin(x)/cos(x)
    assert op2(x, evaluate=False).doit() == sin(x)/cos(x)
