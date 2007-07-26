import py
import sympy as g
from sympy import Symbol, Rational, sin, log, exp, O


def testseries1():
    x = Symbol("x")
    e = sin(x)
    assert e.series(x,0) != 0
    assert e.series(x,0) == O(1, x)
    assert e.series(x,1) == O(x, x)
    assert e.series(x,2) == x + O(x**2, x)
    assert e.series(x,3) == x + O(x**3, x)
    assert e.series(x,4) == x-x**3/6 + O(x**4, x)

    e = (exp(x)-1)/x
    assert e.series(x,2) == 1+x/2+O(x**2, x)
    #this tests, that the Basic.series() cannot do it (but Mul.series can)
    #py.test.raises(g.core.power.pole_error, g.core.basic.Basic.series, e,x,0)

    assert x.series(x,0) == O(1, x)
    assert x.series(x,1) == O(x, x)
    assert x.series(x,2) == x

def testseriesbug1():
    x = Symbol("x")
    assert (1/x).series(x,3) == 1/x
    assert (x+1/x).series(x,3) == x+1/x

def testseries2():
    x = Symbol("x")
    assert ((x+1)**(-2)).series(x,4) == 1-2*x+3*x**2-4*x**3+O(x**4, x)
    assert ((x+1)**(-1)).series(x,4) == 1-x+x**2-x**3+O(x**4, x)
    assert ((x+1)**0).series(x,3) == 1
    assert ((x+1)**1).series(x,3) == 1+x
    assert ((x+1)**2).series(x,3) == 1+2*x+x**2+O(x**3, x)
    assert ((x+1)**3).series(x,3) == 1+3*x+3*x**2+O(x**3, x)

    assert (1/(1+x)).series(x,4) == 1-x+x**2-x**3+O(x**4, x)
    assert (x+3/(1+2*x)).series(x,4) == 3-5*x+12*x**2-24*x**3+O(x**4, x)

    assert ((1/x+1)**3).series(x,3)== 1+x**(-3)+3*x**(-2)+3/x+O(x**3, x)
    assert (1/(1+1/x)).series(x,4) == x-x**2+x**3-O(x**4, x)
    assert (1/(1+1/x**2)).series(x,6) == x**2-x**4+O(x**6, x)

def _testfind():
    # XXX find() doesn't exist
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    p = Rational(5)
    e = a*b+b**p
    assert e.find(b)
    assert not e.find(c)

def xtest_log():
    "too difficult"
    x = Symbol("x")
    ec = exp(Rational(1))
    e = (log(1/x+ec)-ec) / (x*log(1/x+1))
    print e.eval()
    d = e.diff(x)

def _test_bug2(): ### 1/log(0) * log(0) problem
    w = Symbol("w")
    e = (w**(-1)+w**(-log(3)*log(2)**(-1)))**(-1)*(3*w**(-log(3)*log(2)**(-1))+2*w**(-1))
    e = e.eval().expand()
    #should be 3, but is 2
    assert e.series(w,4).subs(w,0)==3

def test_exp():
    x = Symbol("x")
    e = (1+x)**(1/x)
    assert e.series(x,2) == exp(1) - x*exp(1)/2 + O(x**2, x)

def _test_exp2():
    # XXX Fails only when py.test is executed at top level, but passes
    #     when executed within directory (probably caching error)
    x = Symbol("x")
    w = Symbol("w")
    e = w**(1-log(x)/(log(2) + log(x)))
    assert e.eval().series(w,1) != 0

def test_bug3():
    x = Symbol("x")
    e = (2/x+3/x**2)/(1/x+1/x**2)
    assert e.series(x,1) == 3+O(x)

def test_generalexponent():
    x = Symbol("x")
    p = 2
    e = (2/x+3/x**p)/(1/x+1/x**p)
    assert e.eval().series(x,1).leadterm(x) == (3,0)
    p = Rational(1,2)
    e = (2/x+3/x**p)/(1/x+1/x**p)
    assert e.eval().series(x,1).leadterm(x) == (2,0)
    p = Rational(3,2)
    e = (2/x+3/x**p)/(1/x+1/x**p)
    assert e.eval().series(x,1).leadterm(x) == (3,0)

    e=1+x**Rational(1,2)
    assert e.eval().series(x,4) == 1+x**Rational(1,2)

    e=1/(1+x**Rational(1,2))
    assert e.eval().series(x,2) == \
                1+x-x**Rational(1,2)-x**Rational(3,2)+O(x**2, x)

def _test_subsbug1():
    # XXX pole_error doesn't exist
    x = Symbol("x")
    e = 1 + x**Rational(1,2)
    e = e.diff(x)
    py.test.raises(g.core.power.pole_error, e.subs, x, Rational(0))

def test_seriesbug2():
    w = Symbol("w")
    #simple case (1):
    e = ((2*w)/w)**(1+w)
    assert e.series(w,1) == 2 + O(w, w)
    assert e.series(w,1).subs(w,0) == 2

    #test sin
    e = sin(2*w)/w
    assert e.series(w,2) == 2 + O(w**2, w)

    #more complicated case, but sin(x)~x, so the result is the same as in (1)
    e=(sin(2*w)/w)**(1+w)
    assert e.series(w,1) == 2 + O(w)
    assert e.series(w,3) == 2-Rational(4,3)*w**2+w**2*log(2)**2+2*w*log(2)+O(w**3, w)
    assert e.series(w,2).subs(w,0) == 2

def _test_seriesbug3():
    # XXX Fails with infinte recursion
    x = Symbol("x")
    w = Symbol("w")

    #some limits need this series expansion to work:
    e = (w**(-log(5)/log(3))-1/w)**(1/x)
    assert e.series(w,1).subs(log(w),-log(3)*x).subs(w,0) == 5

def test_order():
    #these tests basically define the O's behavior. So be very careful
    #when changing them (as with every test in SymPy :)
    x = Symbol("x")
    assert O(x) == O(x)
    assert O(x**2) == O(x**2)
    assert O(x) != O(x**2)

    assert O(x) + O(x) == O(x)
    assert O(x) - O(x) == O(x)

    assert O(2*x) == O(x)
    assert 2*O(x) == O(x)
    assert O(3*x) == O(x)
    assert 3*O(x) == O(x)
    assert O(x) == O(x*8)
    assert O(x) == O(x)*8

    assert O(x+1) == O(1,x)
    assert O(x)+0 == O(x)
    assert O(x)+0 == 0+O(x)
    assert O(x)+1 != O(x)
    assert O(x)+1 == 1+O(x)
    assert O(x)+x == O(x)
    assert O(x)+x**2 == O(x)
    assert O(x)+x**3 == O(x)

    assert x*O(x) == O(x**2)
    assert O(x)*x == O(x**2)
    assert O(x**3)*x == O(x**4)

    assert ((x+O(x)) - (x+O(x))) == O(x)
    assert O(x)/x == O(1,x)

    assert O(x)*O(x) == O(x**2)

    # XXX these three fail
    #assert O(1,x).diff(x) == O(1,x)
    #assert O(x).diff(x) == O(1,x)
    #assert O(x**2).diff(x) == O(x)

    assert O(x)*Symbol("m") == O(x)
    a = Rational(1,3)+x**(-2)+O(x)
    b = Rational(1,6)-x**(-2)+O(x)
    assert a+b == Rational(1,2) + O(x)

    x = Symbol("w")
    assert O(x)+1 != O(x)

    assert (2+O(x)) != 2
    # XXX These two fail because removeO is not a method 
    #assert (2+O(x)).removeO() == 2
    #assert (2+x+O(x**2)).removeO() == x+2

def test_order_bug():
    x = Symbol("x")
    a = -4
    b = -3/x
    e1 = O(x)*a+O(x)*b
    e2 = O(x)*(a+b)
    assert e1 == O(1,x)
    assert e1 == e2

    assert O(x**2)*(1+2/x+3/x**2) == O(1,x)
    assert O(1+2/x+3/x**2) == O(1/x**2)

def test_order_expand_bug():
    x = Symbol("x")

    assert O(x**2)+x+O(x) == O(x)
    a = O(x**2) + 2*x
    b = 3+O(x)
    assert a+b == 3+O(x)

    e = (2/x+3*x**(-2))*(O(x**3)+x**2)
    assert e.expand() == 3 + O(x)

def _test_expbug4():
    # XXX O(x).series not implemented
    x = Symbol("x")
    assert (log(sin(2*x)/x)*(1+x)).series(x,2) == log(2) + x*log(2) + O(x**2, x)
    assert exp(log(2)+O(x)).series(x,2) == 2 +O(x**2, x)
    assert exp(log(sin(2*x)/x)*(1+x)).series(x,2) == 2 + O(x**2, x)
    assert ((2+O(x))**(1+x)).series(x,2) == 2 + O(x**2, x)

def _test_logbug4():
    # XXX O(x).series not implemented
    x = Symbol("x")
    assert log(2+O(x)).series(x,2) == log(2) + O(x, x)

def _test_expbug5():
    # XXX O(x).series not implemented
    x = Symbol("x")
    assert exp(O(x)).series(x,2) == 1 + O(x**2, x)
    assert exp(log(1+x)/x).series(x,2) == exp(1) + O(x, x)

def test_sinsinbug():
    x = Symbol("x")
    assert sin(sin(x)).series(x,8) == x-x**3/3+x**5/10-8*x**7/315+O(x**8)

def test_issue159():
    x = Symbol("x")
    a=x/(exp(x)-1)
    assert a.series(x,5) == 1 - x/2 - x**4/720 + x**2/12 + O(x**5)

def test_issue105():
    x = Symbol("x")
    f = sin(x**3)**Rational(1,3)
    assert f.series(x, 17) == x - x**7/18 - x**13/3240 + O(x**17)
