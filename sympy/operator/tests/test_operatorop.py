from sympy import Add, cos, Function, sin, symbols
from sympy.abc import x,y,z
from sympy.operator import Operator

f,g,h = symbols('f, g, h', cls=Function)

def test_OperatorAdd():
    assert (1+f)(x) == 1+f(x)
    assert (1+sin+f+g(x))(x) == 1+sin(x)+f(x)+g(x)
    assert (sin-cos)(x) == sin(x)-cos(x)
    assert (sin+1-sin)(x) == 1

def test_OperatorMul():
    assert (1*f)(x) == f(x)
    assert (2*cos)(x) == 2*cos(x)
    assert (2*sin*f*g(x))(x) == 2*sin(x)*f(x)*g(x)
    assert (sin*cos)(x) == sin(x)*cos(x)

def test_OperatorPow():
    assert (sin**2)(x) == sin(x)**2
    assert (sin/cos)(x) == sin(x)/cos(x)
    assert (sin**cos)(x) == sin(x)**cos(x)
    assert ((sin*2)**2)(x) == 4*(sin(x)**2)
    assert ((sin*cos)**2)(x) == (sin(x)*cos(x))**2
    assert ((sin+cos)**2)(x) == (sin(x)+cos(x))**2
    assert ((sin**cos)**2)(x) == (sin(x)**cos(x))**2
