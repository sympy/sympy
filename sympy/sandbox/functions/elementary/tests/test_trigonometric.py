from sympy.sandbox import *

def test_sin():
    x, y = Symbol('x'), Symbol('y')
    assert sin(0) == 0
    assert sin(1) == sin(1)
    assert sin(-1) == -sin(1)
    assert sin(x) == sin(x)
    assert sin(-x) == -sin(x)
    assert sin(pi) == 0
    assert sin(-pi) == 0
    assert sin(2*pi) == 0
    assert sin(-2*pi) == 0
    assert sin(-3*10**73*pi) == 0
    assert sin(7*10**103*pi) == 0
    assert sin(pi/2) == 1
    assert sin(-pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(7*pi/2) == -1
    assert sin(pi/3) == sqrt(3)/2
    assert sin(-2*pi/3) == -sqrt(3)/2
    assert sin(pi/4) == sqrt(2)/2
    assert sin(-pi/4) == -sqrt(2)/2
    assert sin(17*pi/4) == sqrt(2)/2
    assert sin(-3*pi/4) == -sqrt(2)/2
    half = Rational(1,2)
    assert sin(pi/6) == half
    assert sin(-pi/6) == -half
    assert sin(7*pi/6) == -half
    assert sin(-5*pi/6) == -half
    assert sin(pi/105) == sin(pi/105)
    assert sin(-pi/105) == -sin(pi/105)
    assert sin(2 + 3*I) == sin(2 + 3*I)
