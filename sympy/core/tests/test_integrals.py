import py

from sympy import Integral, Symbol, Rational, cos, sin, log, exp
#from sympy.modules.integrals import IntegralError

x = Symbol("x")
t = Symbol("t")
a = Symbol("a")


def integrate(*args, **kwargs):
    ret = expr = args[0]
    for item in args[1:]:
        if isinstance(item, tuple):
            ret = Integral(ret, item[0]==[item[1], item[2]])
        else:
            ret = Integral(ret, item)
    return ret

def test_basics():
    e=(t+1)**2
    assert integrate(e, (t,0,x)).diff(x)==(1+x)**2
    assert integrate(e, (t,0,x)).diff(a)==0

    # XXX Things are automatically evaluated
    #py.test.raises(IntegralError,"integrate(e,(t,0,x), evaluate=False).diff(t)")

    assert integrate(e, (t,a,x)).diff(x)==(1+x)**2
    assert integrate(e, (t,a,x)).diff(x)!=-(1+x)**2
    assert integrate(e, (t,x,a)).diff(x)==-(1+x)**2

    assert integrate(t**2, (t,x,2*x)).diff(x)==7*x**2


def test_integration():
    assert integrate(0, (t,0,x)) == 0
    assert integrate(3, (t,0,x)) == 3*x
    assert integrate(t, (t,0,x)) == x**2/2
    assert integrate(3*t, (t,0,x))== 3*x**2/2
    assert integrate(3*t**2, (t,0,x)) == x**3
    assert integrate(-1/t**2, (t,1,x)) == 1/x-1
    #assert integrate(1/t, (t,1,x)) == log(abs(x))
    assert integrate(t**2+5*t-8, (t,0,x))==x**3/3+5*x**2/2-8*x
    assert integrate(x**2, x) == x**3/3

    b = Symbol("b")
    c = Symbol("c")
    assert integrate(a*t, (t,0,x))==a*x**2/2
    assert integrate(a*t**4, (t,0,x))==a*x**5/5
    assert integrate(a*t**2+b*t+c, (t,0,x))==a*x**3/3+b*x**2/2+c*x

def test_multiple_integration():
    y = Symbol('y')
    assert integrate((x**2)*(y**2), (x,0,1), (y,-1,2)) == Rational(1)
    assert integrate((y**2)*(x**2), x, y) == Rational(1,9)*(x**3)*(y**3)

def _test_integration_table():
    # XXX Fails, gets results oo
    assert integrate(1/(x+1), x) == log(abs(x+1))
    assert integrate(-4*sin(4*x), x) == cos(4*x)
    assert integrate(3*cos(4*x), x) == 3*sin(4*x)/4
    assert integrate(log(x), x) == x*log(x) - x
    assert integrate(log(3*x), x) == x*log(3*x) - x

    # This test is to ensure that the integral table does not look
    # up the *wrong* answer.
    # XXX IntegralError doesn't exist
    #py.test.raises(IntegralError, "integrate(exp(x**x) / (x+1), x)")

    # XXX Fails, doesn't integrate
    assert integrate(exp(x), x) == exp(x)
    assert integrate(2*exp(3*x), x) == 2*exp(3*x)/3

    # XXX No simplify yet
    #from sympy import simplify, diff
    #assert simplify(diff(integrate(x**3 * exp(4*x), x), x)) == x**3 * exp(4*x)
