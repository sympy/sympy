from sympy import (limit, exp, oo, log, sqrt, Limit, sin, floor, cos, ceiling,
                   atan, gamma, Symbol, S, pi, Integral, cot, Rational, I, zoo,
                   tan, cot, integrate, Sum, sign)

from sympy.abc import x, y, z
from sympy.utilities.pytest import XFAIL, raises
from sympy.utilities.iterables import cartes

def test_basic1():
    assert limit(x, x, oo) == oo
    assert limit(x, x, -oo) == -oo
    assert limit(-x, x, oo) == -oo
    assert limit(x**2, x, -oo) == oo
    assert limit(-x**2, x, oo) == -oo
    assert limit(x*log(x), x, 0, dir="+") == 0
    assert limit(1/x, x, oo) == 0
    assert limit(exp(x), x, oo) == oo
    assert limit(-exp(x), x, oo) == -oo
    assert limit(exp(x)/x, x, oo) == oo
    assert limit(1/x - exp(-x), x, oo) == 0
    assert limit(x + 1/x, x, oo) == oo
    assert limit(x - x**2, x, oo) == -oo
    assert limit((1 + x)**(1 + sqrt(2)),x,0) == 1
    assert limit((1 + x)**oo, x, 0) == oo
    assert limit((1 + x)**oo, x, 0, dir='-') == 0
    assert limit((1 + x + y)**oo, x, 0, dir='-') == (1 + y)**(oo)
    assert limit(y/x/log(x), x, 0) == -y*oo
    assert limit(cos(x + y)/x, x, 0) == sign(cos(y))*oo
    raises(NotImplementedError, 'limit(Sum(1/x, (x, 1, y)) - log(y), y, oo)')
    assert limit(Sum(1/x, (x, 1, y)) - 1/y, y, oo) == Sum(1/x, (x, 1, oo))
    assert limit(gamma(1/x + 3), x, oo) == 2

    # approaching 0
    # from dir="+"
    assert limit(1 + 1/x, x, 0) == oo
    # from dir='-'
    # Add
    assert limit(1 + 1/x, x, 0, dir='-') == -oo
    # Pow
    assert limit(x**(-2), x, 0, dir='-') == oo
    assert limit(x**(-3), x, 0, dir='-') == -oo
    assert limit(x**(-Rational(1, 2)), x, 0, dir='-') == (-oo)*I
    assert limit(x**2, x, 0, dir='-') == 0
    assert limit(x**(Rational(1, 2)), x, 0, dir='-') == 0
    assert limit(x**-pi, x, 0, dir='-') == zoo
    assert limit((1 + cos(x))**oo, x, 0) == oo

def test_basic2():
    assert limit(x**x, x, 0, dir="+") == 1
    assert limit((exp(x)-1)/x, x, 0) == 1
    assert limit(1 + 1/x, x, oo) == 1
    assert limit(-exp(1/x), x, oo) == -1
    assert limit(x + exp(-x), x, oo) == oo
    assert limit(x + exp(-x**2), x, oo) == oo
    assert limit(x + exp(-exp(x)), x, oo) == oo
    assert limit(13 + 1/x - exp(-x), x, oo) == 13

def test_basic3():
    assert limit(1/x, x, 0, dir="+") == oo
    assert limit(1/x, x, 0, dir="-") == -oo

def test_basic4():
    assert limit(2*x + y*x, x, 0) == 0
    assert limit(2*x + y*x, x, 1) == 2+y
    assert limit(2*x**8 + y*x**(-3), x, -2) == 512 - y/8
    assert limit(sqrt(x + 1) - sqrt(x), x, oo)==0
    assert integrate(1/(x**3+1),(x,0,oo)) == 2*pi*sqrt(3)/9

def test_issue786():
    assert limit(x*y + x*z, z, 2) == x*y+2*x

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

def test_issue2065():
    assert limit(x**0.5, x, oo) == oo**0.5 == oo
    assert limit(x**0.5, x, 16) == S(16)**0.5
    assert limit(x**0.5, x, 0) == 0
    assert limit(x**(-0.5), x, oo) == 0
    assert limit(x**(-0.5), x, 4) == S(4)**(-0.5)

def test_issue2084():
    # using list(...) so py.test can recalculate values
    tests = list(cartes([x, -x],
                        [-1, 1],
                        [2, 3, Rational(1, 2), Rational(2, 3)],
                        ['-', '+']))
    results = (oo, oo, -oo, oo, -oo*I, oo, -oo*(-1)**Rational(1, 3), oo,
               0, 0, 0, 0, 0, 0, 0, 0,
               oo, oo, oo, -oo, oo, -oo*I, oo, -oo*(-1)**Rational(1, 3),
               0, 0, 0, 0, 0, 0, 0, 0)
    assert len(tests) == len(results)
    for i, (args, res) in enumerate(zip(tests, results)):
        y, s, e, d = args
        eq=y**(s*e)
        try:
            assert limit(eq, x, 0, dir=d) == res
        except AssertionError:
            if 0: # change to 1 if you want to see the failing tests
                print
                print i, res, eq, d, limit(eq, x, 0, dir=d)
            else:
                assert None

def test_issue2085():
    assert limit(sin(x)/x, x, oo) == 0
    assert limit(atan(x), x, oo) == pi/2
    assert limit(gamma(x), x, oo) == oo
    assert limit(cos(x)/x, x, oo) == 0
    assert limit(gamma(x), x, Rational(1, 2)) == sqrt(pi)

@XFAIL
def test_issue2130():
    assert limit((1+y)**(1/y) - S.Exp1, y, 0) == 0

def test_issue1447():
    # using list(...) so py.test can recalculate values
    from sympy import sign
    tests = list(cartes([cot, tan],
                        [-pi/2, 0, pi/2, pi, 3*pi/2],
                        ['-', '+']))
    results = (0, 0, -oo, oo, 0, 0, -oo, oo, 0, 0,
               oo, -oo, 0, 0, oo, -oo, 0, 0, oo, -oo)
    assert len(tests) == len(results)
    for i, (args, res) in enumerate(zip(tests, results)):
        f, l, d= args
        eq=f(x)
        try:
            assert limit(eq, x, l, dir=d) == res
        except AssertionError:
            if 0: # change to 1 if you want to see the failing tests
                print
                print i, res, eq, l, d, limit(eq, x, l, dir=d)
            else:
                assert None

def test_issue835():
    assert limit((1 + x**log(3))**(1/x), x, 0) == 1
    assert limit((5**(1/x) + 3**(1/x))**x, x, 0) == 5

def test_newissue():
    assert limit(exp(1/sin(x))/exp(cot(x)), x, 0) == 1

def test_extended_real_line():
    assert limit(x - oo, x, oo) == -oo
    assert limit(oo - x, x, -oo) == oo
    assert limit(x**2/(x-5) - oo, x, oo) == -oo
    assert limit(1/(x+sin(x)) - oo, x, 0) == -oo
    assert limit(x - oo + 1/x, x, oo) == -oo
    assert limit(x - oo + 1/x, x, 0) == -oo
    assert limit(oo/x, x, oo) == oo

@XFAIL
def test_order_oo():
    from sympy import C
    x = Symbol('x', positive=True, bounded=True)
    assert C.Order(x)*oo != C.Order(1, x)
    assert limit(oo/(x**2 - 4), x, oo) == oo

def test_issue2337():
    raises(NotImplementedError, 'limit(exp(x*y), x, oo)')
    raises(NotImplementedError, 'limit(exp(-x*y), x, oo)')
