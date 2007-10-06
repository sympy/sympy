
from sympy import *

def test_exp():

    x, y = symbols('xy')

    assert exp(nan) == nan

    assert exp(oo) == oo
    assert exp(-oo) == 0

    assert exp(0) == 1
    assert exp(1) == E

    assert exp(pi*I/2) == I
    assert exp(pi*I) == -1
    assert exp(3*pi*I/2) == -I
    assert exp(2*pi*I) == 1

    assert exp(log(x)) == x
    assert exp(2*log(x)) == x**2
    assert exp(pi*log(x)) == x**pi

    assert exp(17*log(x) + E*log(y)) == x**17 * y**E

    assert exp(x*log(x)) != x**x


def test_log():

    assert log(nan) == nan

    assert log(oo) == oo
    assert log(-oo) == oo

    assert log(0) == -oo

    assert log(1) == 0
    assert log(-1) == I*pi

    assert log(E) == 1
    assert log(-E).expand() == 1 + I*pi

    assert log(pi) == log(pi)
    assert log(-pi).expand() == log(pi) + I*pi

    assert log(17) == log(17)
    assert log(-17) == log(17) + I*pi

    assert log(I) == I*pi/2
    assert log(-I) == -I*pi/2

    assert log(17*I) == I*pi/2 + log(17)
    assert log(-17*I).expand() == -I*pi/2 + log(17)

    assert log(oo*I) == oo
    assert log(-oo*I) == oo

    assert exp(-log(3))**(-1) == 3

    x, y = symbols('xy')

    assert log(x) == log(x)
    assert log(x*y) != log(x) + log(y)

    #assert log(x**2) != 2*log(x)
    assert log(x**2).expand() == 2*log(x)
    assert log(x**y) != y*log(x)

    #I commented this test out, because it doesn't work well with caching and
    #thus completely breaks limits, that rely on log(exp(x)) -> x
    #simplification. --Ondrej
    #assert log(exp(x)) != x

    x, y = symbols('xy', real=True)

    assert log(x) == log(x)
    #assert log(x*y) != log(x) + log(y)
    assert log(x*y).expand() == log(x) + log(y)

    #assert log(x**2) != 2*log(x)
    assert log(x**2).expand() == 2*log(x)
    assert log(x**y) != y*log(x)

    assert log(exp(x)) == x
    #assert log(-exp(x)) != x + I*pi
    assert log(-exp(x)).expand() == x + I*pi

    k = Symbol('k', positive=True)

    assert log(-x) == log(-x)
    assert log(-k) == log(-k)

    assert log(x, 2) == log(x)/log(2)
    assert log(E, 2) == 1/log(2)

def test_log_apply_evalf():
    value = (log(3)/log(2)-1).evalf()

    assert value.epsilon_eq(Real("0.58496250072115618145373"))

def test_log_simplify():
    x = Symbol("x")
    assert log(x**2) == 2*log(x)
    assert log(x**(2+log(2))) == (2+log(2))*log(x)
