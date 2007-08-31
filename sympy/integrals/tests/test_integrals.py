
from sympy import *

x = Symbol('x')

def test_proper_integral():
    pass

def test_improper_integral():
    assert integrate(ln(x), (x,0,1)) == -1
    assert integrate(x**(-2), (x,1,oo)) == 1

'''
import py
from sympy import *
from sympy.integrals import integrate, IntegralError

x = Symbol("x")
t = Symbol("t")
a = Symbol("a")

def diff(expr, sym):
    from sympy.core import Derivative
    if isinstance(sym, Symbol):
        return Derivative(expr, sym)
    else:
        return Derivative(expr, sym[0]==[sym[1], sym[2]])


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
    assert integrate(-1/t**2, (t,1,x)) == 1/x-1
    assert integrate(1/t, (t,1,x)) == log(abs(x))
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

def test_integration_table():
    assert integrate(log(x), x) == x*log(x) - x
    assert integrate(log(3*x), x) == x*log(3*x) - x
    assert integrate(log(x+1), x) == (x+1)*log(x+1)-x
    assert integrate(1/(x+1), x) == log(abs(x+1))
    assert integrate(5*x**5/(2*x**6+5), x) == 5*log(abs(2*x**6+5))/12

    # XXX Equality testing fails
    #assert integrate((3*x+2)**a, x) == (3*x+2)**(a+1) / (3*a+3)
    assert integrate((3*x+2)**a, x) == (3*x+2)**(a+1) / (1+a) / 3
    assert integrate(x/(2*x+1), x) == (2*x-log(abs(2*x+1)))/4
    assert integrate(2*x/(4*x+a)**2, x) == a/(a+4*x)/8 + log(abs(a+4*x))/8

    assert integrate(-4*sin(4*x), x) == cos(4*x)
    assert integrate(3*cos(4*x), x) == 3*sin(4*x)/4
    assert integrate(3*t*cos(4*x) + 4*sin(x), x) == 3*t*sin(4*x)/4 - 4*cos(x)

    # This test is to ensure that the integral table does not look
    # up the *wrong* answer.
    py.test.raises(IntegralError, "integrate(exp(x**x) / (x+1), x)")

    assert integrate(exp(x), x) == exp(x)
    assert integrate(2*exp(3*x), x) == 2*exp(3*x)/3
    assert diff(integrate(x**3 * exp(4*x), x), x) == x**3 * exp(4*x)

def test_improper_integral():
    assert integrate(ln(x), (x,0,1)) == -1
    assert integrate(x**(-2), (x,1,oo)) == 1
'''