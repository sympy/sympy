from sympy import symbols, integrate, Integral, exp, oo, Symbol, Rational, \
    log, sin, cos, pi, E, I, Polynomial, LambertW, diff
from sympy.utilities.pytest import XFAIL
from sympy.physics.units import m, s
import py
from py.test import skip

x,y,a,t = symbols('xyat')
n = Symbol('n', integer=True)

def test_improper_integral():
    assert integrate(log(x), (x, 0, 1)) == -1
    assert integrate(x**(-2), (x, 1, oo)) == 1

def test_basics():
    e=(t+1)**2
    assert integrate(e, (t,0,x), evaluate=False).diff(x)==((1+x)**2).expand()
    assert integrate(e, (t,0,x), evaluate=False).diff(a)==0

    assert integrate(e, (t,a,x), evaluate=False).diff(x)==((1+x)**2).expand()
    assert integrate(e, (t,x,a), evaluate=False).diff(x)==(-(1+x)**2).expand()

    assert integrate(t**2, (t,x,2*x), evaluate=False).diff(x)==7*x**2

@XFAIL
def test_unevaluated():
    py.test.raises(IntegralError,"integrate(e, (t,0,x), evaluate=False).diff(t)")
    py.test.raises(IntegralError,"integrate(e, t, evaluate=False).diff(x)")



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
    assert integrate((x**2)*(y**2), (x,0,1), (y,-1,2)) == Rational(1)
    assert integrate((y**2)*(x**2), x, y) == Rational(1,9)*(x**3)*(y**3)
    assert integrate(1/(x+3)/(1+x)**3, x) == \
        -(1 + x)**(-2)/4 - log(3 + x)/8 + 1/(1 + x)/4 + log(1 + x)/8

def test_issue433():
    assert integrate(exp(-x), (x,0,oo)) == 1

def test_issue461():
    assert integrate(x**Rational(3,2), x) == 2*x**Rational(5,2)/5
    assert integrate(x**Rational(1,2), x) == 2*x**Rational(3,2)/3
    assert integrate(x**Rational(-3,2), x) == -2*x**Rational(-1,2)

def test_integrate_poly():
    p = Polynomial(x + x**2*y + y**3, var=[x,y])
    qx= integrate(p, x)
    qy= integrate(p, y)

    assert isinstance(qx, Polynomial) == True
    assert isinstance(qy, Polynomial) == True

    assert qx.var == [x,y]
    assert qy.var == [x,y]

    assert qx == x**2/2 + x**3*y/3 + x*y**3
    assert qy == x*y + x**2*y**2/2 + y**4/4


@XFAIL #depends on poly.subs -> poly
def test_integrate_poly_defined():
    p = Polynomial(x + x**2*y + y**3, var=[x,y])
    Qx= integrate(p, (x,0,1))
    Qy= integrate(p, (y,0,pi))

    assert isinstance(Qx, Polynomial) == True
    assert isinstance(Qy, Polynomial) == True

    assert Qx.var == [y]
    assert Qy.var == [x]

    assert Qx == Rational(1,2) + y/3 + y**3
    assert Qy == pi**4/4 + pi*x + pi**2*x**2/2


def test_integrate_varommited():
    y = Symbol('y')
    assert integrate(2)     == 2
    assert integrate(x)     == x**2/2
    assert integrate(x*y)   == x**2*y**2/4


def test_integrate_poly_accurately():
    y = Symbol('y')
    assert integrate(x*sin(y), x)       == x**2*sin(y)/2

    # when passed to risch_norman, this will be a CPU hog, so this really
    # checks, that integrated function is recognized as polynomial
    assert integrate(x**1000*sin(y), x) == x**1001*sin(y)/1001

def test_issue536():
    y = Symbol('y')
    assert integrate(x**2, y) == x**2*y
    assert integrate(x**2, (y, -1, 1)) == 2*x**2

# works in sympy and py.test but hangs in `setup.py test`
def test_integrate_linearterm_pow():
    # check integrate((a*x+b)^c, x)  --  #400
    y = Symbol('y')
    assert integrate(x**y, x) == x**(y+1)/(y+1)
    assert integrate((exp(y)*x + 1/y)**(1+sin(y)), x)   == exp(-y)*(exp(y)*x + 1/y)**(2+sin(y)) / (2+sin(y))

def test_issue519():
    assert integrate(pi*x**Rational(1,2),x) == 2*pi*x**Rational(3,2)/3
    assert integrate(pi*x**Rational(1,2) + E*x**Rational(3,2),x) == \
                                               2*pi*x**Rational(3,2)/3  + \
                                               2*E *x**Rational(5,2)/5
def test_issue524():
    assert integrate(cos((n+1) * x), x)   == sin(x*(n+1)) / (n+1)
    assert integrate(cos((n-1) * x), x)   == sin(x*(n-1)) / (n-1)

    assert integrate(cos((n+1) * x) + cos((n-1) * x), x) == \
                                             sin(x*(n+1)) / (n+1)  + \
                                             sin(x*(n-1)) / (n-1)

def test_issue565():
    assert integrate(-1./2 * x * sin(n * pi * x/2), [x, -2, 0])  == 2*cos(pi*n)/(pi*n)
    assert integrate(-Rational(1)/2 * x * sin(n * pi * x/2), [x, -2, 0]) \
                                                                 == 2*cos(pi*n)/(pi*n)


def test_rational_functions():
    half = Rational(1,2)
    integrate(1/(x**2+x+1), x) == -I*3**half*log(half + x - half*I*3**half)/3 +\
                                   I*3**half*log(half + x + half*I*3**half)/3

    integrate(1/(x**3+1), x) == log(1 + x)/3 - \
        (Rational(1,6) - I*3**half/6)*log(half - x - I*3**half/2) - \
        (Rational(1,6) + I*3**half/6)*log(half - x + I*3**half/2)

def test_issue587(): # remove this when fresnel itegrals are implemented
    assert integrate(sin(x**2), x) == Integral(sin(x**2), x)

def test_integrate_units():
    assert integrate(x * m/s, (x, 1*s, 5*s)) == 12*m*s

def test_transcendental_functions():
    assert integrate(LambertW(2*x), x) == -x + x*LambertW(2*x) + x/LambertW(2*x)

def test_issue641():
    f=4*log(x)-2*log(x)**2
    fid=diff(integrate(f,x),x)
    assert abs(f.subs(x,42).evalf() - fid.subs(x,42).evalf()) < 1e-10
