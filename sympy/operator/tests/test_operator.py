from sympy import cos, Function, sin, symbols
from sympy.abc import x,y,z
from sympy.operator import Op

f,g,h = symbols('f, g, h', cls=Function)

def test_Op():
    assert Op(1)(x,y,z) == 1
    assert Op(sin)(x) == sin(x)
    assert Op(cos)(y) == cos(y)
    assert Op(f)(x) == f(x)
    assert Op(g)(x,y) == g(x,y)
    assert Op(h(x))(x,y) == h(x)

    assert Op(sin)(x,y) == sin(x)
    assert Op(sin, (1,))(x,y,z) == sin(y)

def test_AppliedOp1():
    """
    Cases where AppliedOp is force-evaluated.
    """
    assert Op(1)(x,y, evaluate=False) == 1
    assert Op(f(x))(y,z, evaluate=False) == f(x)
    assert Op(sin)(x, evaluate=False) == sin(x)

def test_AppliedOp2():
    """
    Cases where AppliedOp is not force-evaluated.
    """
    op1 = sin+cos
    op2 = sin/cos
    assert op1(x, evaluate=False) != sin(x)+cos(x)
    assert op1(x, evaluate=False).doit() == sin(x)+cos(x)
    assert op2(x, evaluate=False) != sin(x)/cos(x)
    assert op2(x, evaluate=False).doit() == sin(x)/cos(x)
