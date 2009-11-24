from sympy import limit, exp, oo, log, sqrt, Limit, sin, floor, cos, ceiling, \
        atan, Symbol, S, pi, Integral, cot
from sympy.abc import x, y, z
from sympy.utilities.pytest import XFAIL

def test_basic1():
    assert limit(x, x, oo) == oo
    assert limit(x, x, -oo) == -oo
    assert limit(-x, x, oo) == -oo
    assert limit(x**2, x, -oo) == oo
    assert limit(-x**2, x, oo) == -oo
    assert limit(x*log(x), x, 0, dir="+") == 0
    assert limit(1/x,x,oo) == 0
    assert limit(exp(x),x,oo) == oo
    assert limit(-exp(x),x,oo) == -oo
    assert limit(exp(x)/x,x,oo) == oo
    assert limit(1/x-exp(-x),x,oo) == 0
    assert limit(x+1/x,x,oo) == oo
    assert limit(x-x**2,x,oo) == -oo


def test_basic2():
    assert limit(x**x, x, 0, dir="+") == 1
    assert limit((exp(x)-1)/x, x, 0) == 1
    assert limit(1+1/x,x,oo) == 1
    assert limit(-exp(1/x),x,oo) == -1
    assert limit(x+exp(-x),x,oo) == oo
    assert limit(x+exp(-x**2),x,oo) == oo
    assert limit(x+exp(-exp(x)),x,oo) == oo
    assert limit(13+1/x-exp(-x),x,oo) == 13

def test_basic3():
    assert limit(1/x, x, 0, dir="+") == oo
    assert limit(1/x, x, 0, dir="-") == -oo

def test_basic4():
    assert limit(2*x + y*x, x, 0) == 0
    assert limit(2*x + y*x, x, 1) == 2+y
    assert limit(2*x**8 + y*x**(-3), x, -2) == 512-y/8
    assert limit(sqrt(x+1)-sqrt(x),x,oo)==0

def test_issue786():
    assert limit(x*y+x*z, z, 2) == x*y+2*x

def test_Limit():
    assert Limit(sin(x)/x, x, 0) != 1
    assert Limit(sin(x)/x, x, 0).doit() == 1

def test_floor():
    assert limit(floor(x), x, -2, "+") == -2
    assert limit(floor(x), x, -2, "-") == -3
    assert limit(floor(x), x, -1, "+") == -1
    assert limit(floor(x), x, -1, "-") == -2
    assert limit(floor(x), x, 0, "+") == 0
    assert limit(floor(x), x, 0, "-") == -1
    assert limit(floor(x), x, 1, "+") == 1
    assert limit(floor(x), x, 1, "-") == 0
    assert limit(floor(x), x, 2, "+") == 2
    assert limit(floor(x), x, 2, "-") == 1
    assert limit(floor(x), x, 248, "+") == 248
    assert limit(floor(x), x, 248, "-") == 247

    # note: if any of the tests below fails, just comment it out. General fix
    # needs better assumptions handling.

    # this doesn't work, it requires robust assumptions:
    assert limit(floor(sin(x)), x, 0, "+") == 0
    assert limit(floor(sin(x)), x, 0, "-") == -1
    assert limit(floor(cos(x)), x, 0, "+") == 0
    assert limit(floor(cos(x)), x, 0, "-") == 0

    # this doesn't work, it requires robust assumptions:
    assert limit(floor(5+sin(x)), x, 0, "+") == 5
    #assert limit(floor(5+sin(x)), x, 0, "-") == 4
    #assert limit(floor(5+cos(x)), x, 0, "+") == 5
    #assert limit(floor(5+cos(x)), x, 0, "-") == 5

def test_ceiling():
    assert limit(ceiling(x), x, -2, "+") == -1
    assert limit(ceiling(x), x, -2, "-") == -2
    assert limit(ceiling(x), x, -1, "+") == 0
    assert limit(ceiling(x), x, -1, "-") == -1
    assert limit(ceiling(x), x, 0, "+") == 1
    assert limit(ceiling(x), x, 0, "-") == 0
    assert limit(ceiling(x), x, 1, "+") == 2
    assert limit(ceiling(x), x, 1, "-") == 1
    assert limit(ceiling(x), x, 2, "+") == 3
    assert limit(ceiling(x), x, 2, "-") == 2
    assert limit(ceiling(x), x, 248, "+") == 249
    assert limit(ceiling(x), x, 248, "-") == 248

    # note: if any of the tests below fails, just comment it out. General fix
    # needs better assumptions handling.

    # this doesn't work, it requires robust assumptions:
    #assert limit(ceiling(sin(x)), x, 0, "+") == 1
    assert limit(ceiling(sin(x)), x, 0, "-") == 0
    assert limit(ceiling(cos(x)), x, 0, "+") == 1
    assert limit(ceiling(cos(x)), x, 0, "-") == 1

    # this doesn't work, it requires robust assumptions:
    #assert limit(ceiling(5+sin(x)), x, 0, "+") == 6
    assert limit(ceiling(5+sin(x)), x, 0, "-") == 5
    assert limit(ceiling(5+cos(x)), x, 0, "+") == 6
    assert limit(ceiling(5+cos(x)), x, 0, "-") == 6

def test_atan():
    x = Symbol("x", real=True)
    assert limit(atan(x)*sin(1/x), x, 0) == 0
    assert limit(atan(x) + sqrt(x+1) - sqrt(x), x, oo) == pi/2

def test_abs():
    assert limit(abs(x), x, 0) == 0
    assert limit(abs(sin(x)), x, 0) == 0
    assert limit(abs(cos(x)), x, 0) == 1
    assert limit(abs(sin(x+1)), x, 0) == sin(1)

def test_heuristic():
    x = Symbol("x", real=True)
    assert limit(log(2+sqrt(atan(x))*sqrt(sin(1/x))), x, 0) == log(2)

def test_issue772():
    z = Symbol("z", positive=True)
    f = -1/z*exp(-z*x)
    assert limit(f, x, oo) == 0
    assert f.limit(x, oo) == 0

def test_exponential():
    n = Symbol('n')
    assert limit((1+x/n)**n,n,oo) == exp(x)
    assert limit((1+x/(2*n))**n,n,oo) == exp(x/2)
    assert limit((1+x/(2*n+1))**n,n,oo) == exp(x/2)
    assert limit(((x-1)/(x+1))**x,x,oo) == exp(-2)

@XFAIL
def test_exponential2():
    n = Symbol('n')
    assert limit((1+x/(n+sin(n)))**n,n,oo) == exp(x)

def test_doit():
    f = Integral(2 * x, x)
    l = Limit(f, x, oo)
    assert l.doit() == oo

@XFAIL
def test_doit2():
    f = Integral(2 * x, x)
    l = Limit(f, x, oo)
    # limit() breaks on the contained Integral.
    assert l.doit(deep = False) == l

def test_bug693a():
    assert sin(sin(x+1)+1).limit(x,0) == sin(sin(1)+1)

def test_issue693():
    assert limit( (1-cos(x))/x**2, x, S(1)/2) == 4 - 4*cos(S(1)/2)
    assert limit(sin(sin(x+1)+1), x, 0) == sin(1 + sin(1))
    assert limit(abs(sin(x+1)+1), x, 0) == 1 + sin(1)

def test_issue991():
    assert limit(1/(x+3), x, 2) == S(1)/5
    assert limit(1/(x+pi), x, 2) == S(1)/(2+pi)
    assert limit(log(x)/(x**2+3), x, 2) == log(2)/7
    assert limit(log(x)/(x**2+pi), x, 2) == log(2)/(4+pi)

def test_issue1448():
    assert limit(cot(x),x,0,dir='+') == oo
    assert limit(cot(x),x,pi/2,dir='+') == 0
