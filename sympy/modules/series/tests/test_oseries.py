from sympy import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
zero = Basic.Zero()
o = Basic.One()


def test_simple_1():
    assert x.oseries(x) == 0
    assert x.oseries(1/x) == 0
    assert x.oseries(x**2) == x
    assert y.oseries(x) == y
    assert y.oseries(1/x) == 0
    assert (Rational(3,4)).oseries(x) == Rational(3,4)
    assert (Rational(3,4)).oseries(1/x) == 0

def test_mul_0():
    assert (x*ln(x)).oseries(x**5) == x*ln(x)
    assert (x*ln(x)).oseries(x) == x*ln(x)
    assert (x*ln(x)).oseries(x*ln(x)) == 0
    assert (x*ln(x)).oseries(ln(x)) == 0
    assert (x*ln(x)).oseries(x**2*ln(x)) == x*ln(x)
    assert (x*ln(x)).oseries(1/x) == 0

def test_mul_1():
    assert (x*ln(2+x)).oseries(x**5) == x**2/2-x**3/8+x**4/24+x*log(2)
    assert (x*ln(1+x)).oseries(x**5) == x**2-x**3/2+x**4/3

def test_pow_0():
    assert (x**2).oseries(x**5) == x**2
    assert (x**2).oseries(x) == 0
    assert (1/x).oseries(x) == 1/x
    assert (x**(Rational(2,3))).oseries(x) == (x**(Rational(2,3)))
    assert (x**(Rational(3,2))).oseries(x) == 0

def test_pow_1():
    assert ((1+x)**2).oseries(x**5) == 1+2*x+x**2
    assert ((1+x)**2).oseries(x**2) == 1+2*x
    assert ((1+x)**2).oseries(x) == 1
    assert ((1+x)**2).oseries(1/x) == 0

def test_geometric_1():
    assert (1/(1-x)).oseries(x**5) == 1+x+x**2+x**3+x**4
    assert (x/(1-x)).oseries(x**5) == x+x**2+x**3+x**4
    assert (x**3/(1-x)).oseries(x**5) == x**3+x**4

def test_sqrt_1():
    assert sqrt(1+x).oseries(x**5) == 1+x/2-5*x**4/128-x**2/8+x**3/16

def test_exp_1():
    assert exp(x).oseries(x**5) == 1+x+x**2/2+x**3/6+x**4/24
    assert exp(1/x).oseries(x**5) == exp(1/x)
    assert exp(1/(1+x)).oseries(x**4) == (E*(1-x-13*x**3/6+3*x**2/2)).expand()
    assert exp(2+x).oseries(x**5) == (exp(2)*(1+x+x**2/2+x**3/6+x**4/24)).expand()

def test_exp_sqrt_1():
    assert exp(1+sqrt(x)).oseries(x**2) == (exp(1)*(1+sqrt(x)+x/2+sqrt(x)*x/6)).expand()

def test_power_x_x():
    assert (exp(x*ln(x))).oseries(x**3) == 1+x*log(x)+x**2*log(x)**2/2+x**3*log(x)**3/6
    assert (x**x).oseries(x) == 1+x*ln(x)
