from sympy import sin, cos, exp, E, series, oo, S, Derivative, O, Integral, \
                  Function, log, sqrt, Symbol
from sympy.abc import x, y
from sympy.utilities.pytest import raises

def test_sin():
    e1 = sin(x).series(x, 0)
    e2 = series(sin(x), x, 0)
    assert e1 == e2

def test_cos():
    e1 = cos(x).series(x, 0)
    e2 = series(cos(x), x, 0)
    assert e1 == e2

def test_exp():
    e1 = exp(x).series(x, 0)
    e2 = series(exp(x), x, 0)
    assert e1 == e2

def test_exp2():
    e1 = exp(cos(x)).series(x, 0)
    e2 = series(exp(cos(x)),x,0)
    assert e1 == e2


def test_2124():
    assert series(1, x) == 1
    assert S(0).lseries(x).next() == 0
    assert cos(x).series() == cos(x).series(x)
    raises(ValueError, 'cos(x+y).series()')
    raises(ValueError, 'x.series(dir="")')

    assert (cos(x).series(x, 1).removeO() -
            cos(x + 1).series(x).removeO().subs(x, x - 1)).expand() == 0
    e = cos(x).series(x, 1, n=None)
    assert [e.next() for i in range(2)] == [cos(1), (1 - x)*sin(1)]
    e = cos(x).series(x, 1, n=None, dir='-')
    assert [e.next() for i in range(2)] == [cos(1), (1 - x)*sin(1)]
    assert abs(x).series(x, 1, dir='-') == x
    assert exp(x).series(x, 1, dir="-", n=3) == E - E*(1 - x) + E*(1 - x)**2/2

    D = Derivative
    assert D(x**2 + x**3*y**2, x, 2, y, 1).series(x).doit() == 12*x*y
    assert D(cos(x), x).lseries().next() == D(1, x)
    assert D(exp(x), x).series(n=3) == D(1, x) + D(x, x) + D(x**2/2, x) + O(x**3)

    assert Integral(x, (x, 1, 3),(y, 1, x)).series(x) == -4 + 4*x

    assert (1 + x + O(x**2)).getn() == 2
    assert (1 + x).getn() == None

    assert ((1/sin(x))**oo).series() == oo
    assert ((sin(x))**y).series(x) == 0**y

    raises(NotImplementedError, 'series(Function("f")(x))')

    assert sin(1/x).series(x, oo, n=5) == 1/x - 1/(6*x**3)
    assert abs(x).series(x, oo, n=5, dir='+') == x
    assert abs(x).series(x, -oo, n=5, dir='-') == -x
    assert abs(-x).series(x, oo, n=5, dir='+') == x
    assert abs(-x).series(x, -oo, n=5, dir='-') == -x

    assert exp(x*log(x)).series(n=3) == \
           1 + x*log(x) + x**2*log(x)**2/2 + O(x**3)
    # XXX is this right? If not, fix "ngot > n" handling in expr.
    p = Symbol('p', positive=True)
    assert exp(sqrt(p)**3*log(p)).series(n=3) == \
        1 + p**S('3/2')*log(p) + O(p**3)
