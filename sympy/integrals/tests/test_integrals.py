
from sympy import *
from sympy.utilities.pytest import XFAIL

x,a,t = symbols('xat')

def test_proper_integral():
    pass

def test_improper_integral():
    assert integrate(log(x), (x, 0, 1)) == -1
    assert integrate(x**(-2), (x, 1, oo)) == 1


def diff(expr, sym):
    from sympy.core import Derivative
    if isinstance(sym, Symbol):
        return Derivative(expr, sym)
    else:
        return Derivative(expr, sym[0]==[sym[1], sym[2]])


@XFAIL  # ok, but (1+x)**2 == 1+2*x+x**2  fails
def test_basics():
    e=(t+1)**2
    assert integrate(e, (t,0,x), evaluate=False).diff(x)==(1+x)**2
    assert integrate(e, (t,0,x), evaluate=False).diff(a)==0

    py.test.raises(IntegralError,"integrate(e, (t,0,x), evaluate=False).diff(t)")
    py.test.raises(IntegralError,"integrate(e, t, evaluate=False).diff(x)")

    assert integrate(e, (t,a,x), evaluate=False).diff(x)==(1+x)**2
    assert integrate(e, (t,a,x), evaluate=False).diff(x)!=-(1+x)**2
    assert integrate(e, (t,x,a), evaluate=False).diff(x)==-(1+x)**2

    assert integrate(t**2, (t,x,2*x), evaluate=False).diff(x)==7*x**2


def test_integration():
    assert integrate(0, (t,0,x)) == 0
    assert integrate(3, (t,0,x)) == 3*x
    assert integrate(t, (t,0,x)) == x**2/2
    assert integrate(3*t, (t,0,x))== 3*x**2/2
    assert integrate(3*t**2, (t,0,x)) == x**3
    assert integrate(1/t, (t,1,x)) == log(x)
    assert integrate(-1/t**2, (t,1,x)) == 1/x-1
    assert integrate(t**2+5*t-8, (t,0,x)) == x**3/3+5*x**2/2-8*x
    assert integrate(x**2, x) == x**3/3
    assert integrate((3*t*x)**5, x) == (3*t)**5 * x**6 / 6

    b = Symbol("b")
    c = Symbol("c")
    assert integrate(a*t, (t,0,x))==a*x**2/2
    assert integrate(a*t**4, (t,0,x))==a*x**5/5
    assert integrate(a*t**2+b*t+c, (t,0,x))==a*x**3/3+b*x**2/2+c*x

def test_multiple_integration():
    y = Symbol('y')
    assert integrate((x**2)*(y**2), (x,0,1), (y,-1,2)) == Rational(1)
    assert integrate((y**2)*(x**2), x, y) == Rational(1,9)*(x**3)*(y**3)


def test_issue433():
    assert integrate(exp(-x), (x,0,oo)) == 1
