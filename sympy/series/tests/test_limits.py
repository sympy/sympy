from sympy import limit, exp, oo, log, sqrt, Limit, sin
from sympy.abc import x, y, z

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

def test_basic3():
    assert limit(2*x + y*x, x, 0) == 0
    assert limit(2*x + y*x, x, 1) == 2+y
    assert limit(2*x**8 + y*x**(-3), x, -2) == 512-y/8
    assert limit(sqrt(x+1)-sqrt(x),x,oo)==0

def test_issue786():
    assert limit(x*y+x*z, z, 2) == x*y+2*x

def test_Limit():
    assert Limit(sin(x)/x, x, 0) != 1
    assert str(Limit(1/x, x, 0)) == "Limit(1/x, x, 0, dir='+')" != 1
    assert Limit(sin(x)/x, x, 0).doit() == 1
