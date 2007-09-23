
from sympy import *

def test_sin():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert sin(nan) == nan

    assert sin(oo*I) == oo*I
    assert sin(-oo*I) == -oo*I

    assert sin(0) == 0

    assert sin(1) == sin(1)
    assert sin(-1) == -sin(1)

    assert sin(x) == sin(x)
    assert sin(-x) == -sin(x)

    assert sin(pi*I) == sinh(pi)*I
    assert sin(-pi*I) == -sinh(pi)*I

    assert sin(2**1024 * E) == sin(2**1024 * E)
    assert sin(-2**1024 * E) == -sin(2**1024 * E)

    assert sin(pi) == 0
    assert sin(-pi) == 0
    assert sin(2*pi) == 0
    assert sin(-2*pi) == 0
    assert sin(-3*10**73*pi) == 0
    assert sin(7*10**103*pi) == 0

    assert sin(pi/2) == 1
    assert sin(-pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(7*pi/2) == -1

    assert sin(pi/3) == Basic.Half()*sqrt(3)
    assert sin(-2*pi/3) == -Basic.Half()*sqrt(3)

    assert sin(pi/4) == Basic.Half()*sqrt(2)
    assert sin(-pi/4) == -Basic.Half()*sqrt(2)
    assert sin(17*pi/4) == Basic.Half()*sqrt(2)
    assert sin(-3*pi/4) == -Basic.Half()*sqrt(2)

    assert sin(pi/6) == Basic.Half()
    assert sin(-pi/6) == -Basic.Half()
    assert sin(7*pi/6) == -Basic.Half()
    assert sin(-5*pi/6) == -Basic.Half()

    assert sin(pi/105) == sin(pi/105)
    assert sin(-pi/105) == -sin(pi/105)

    assert sin(2 + 3*I) == sin(2 + 3*I)

    assert sin(x*I) == sinh(x)*I

    assert sin(k*pi) == 0
    assert sin(17*k*pi) == 0

    assert sin(k*pi*I) == sinh(k*pi)*I

def test_cos():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert cos(nan) == nan

    assert cos(oo*I) == oo
    assert cos(-oo*I) == oo

    assert cos(0) == 1

    assert cos(1) == cos(1)
    assert cos(-1) == cos(1)

    assert cos(x) == cos(x)
    assert cos(-x) == cos(x)

    assert cos(pi*I) == cosh(pi)
    assert cos(-pi*I) == cosh(pi)

    assert cos(2**1024 * E) == cos(2**1024 * E)
    assert cos(-2**1024 * E) == cos(2**1024 * E)

    assert cos(pi/2) == 0
    assert cos(-pi/2) == 0
    assert cos(pi/2) == 0
    assert cos(-pi/2) == 0
    assert cos((-3*10**73+1)*pi/2) == 0
    assert cos((7*10**103+1)*pi/2) == 0

    assert cos(pi) == -1
    assert cos(-pi) == -1
    assert cos(5*pi) == -1
    assert cos(8*pi) == 1

    assert cos(pi/3) == Basic.Half()
    assert cos(-2*pi/3) == -Basic.Half()

    assert cos(pi/4) == Basic.Half()*sqrt(2)
    assert cos(-pi/4) == Basic.Half()*sqrt(2)
    assert cos(11*pi/4) == -Basic.Half()*sqrt(2)
    assert cos(-3*pi/4) == -Basic.Half()*sqrt(2)

    assert cos(pi/6) == Basic.Half()*sqrt(3)
    assert cos(-pi/6) == Basic.Half()*sqrt(3)
    assert cos(7*pi/6) == -Basic.Half()*sqrt(3)
    assert cos(-5*pi/6) == -Basic.Half()*sqrt(3)

    assert cos(pi/105) == cos(pi/105)
    assert cos(-pi/105) == cos(pi/105)

    assert cos(2 + 3*I) == cos(2 + 3*I)

    assert cos(x*I) == cosh(x)

    assert cos(k*pi) == cos(k*pi)
    assert cos(17*k*pi) == cos(17*k*pi)

    assert cos(k*pi*I) == cosh(k*pi)

def test_tan():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert tan(nan) == nan

    assert tan(oo*I) == I
    assert tan(-oo*I) == -I

    assert tan(0) == 0

    assert tan(1) == tan(1)
    assert tan(-1) == -tan(1)

    assert tan(x) == tan(x)
    assert tan(-x) == -tan(x)

    assert tan(pi*I) == tanh(pi)*I
    assert tan(-pi*I) == -tanh(pi)*I

    assert tan(2**1024 * E) == tan(2**1024 * E)
    assert tan(-2**1024 * E) == -tan(2**1024 * E)

    assert tan(pi) == 0
    assert tan(-pi) == 0
    assert tan(2*pi) == 0
    assert tan(-2*pi) == 0
    assert tan(-3*10**73*pi) == 0
    assert tan(7*10**103*pi) == 0

    assert tan(pi/2) == tan(pi/2)
    assert tan(-pi/2) == -tan(pi/2)
    assert tan(5*pi/2) == tan(5*pi/2)
    assert tan(7*pi/2) == tan(7*pi/2)

    assert tan(pi/3) == sqrt(3)
    assert tan(-2*pi/3) == sqrt(3)

    assert tan(pi/4) == Basic.One()
    assert tan(-pi/4) == -Basic.One()
    assert tan(17*pi/4) == Basic.One()
    assert tan(-3*pi/4) == Basic.One()

    assert tan(pi/6) == 1/sqrt(3)
    assert tan(-pi/6) == -1/sqrt(3)
    assert tan(7*pi/6) == 1/sqrt(3)
    assert tan(-5*pi/6) == 1/sqrt(3)

    assert tan(pi/105) == tan(pi/105)
    assert tan(-pi/105) == -tan(pi/105)

    assert tan(2 + 3*I) == tan(2 + 3*I)

    assert tan(x*I) == tanh(x)*I

    assert tan(k*pi) == 0
    assert tan(17*k*pi) == 0

    assert tan(k*pi*I) == tanh(k*pi)*I

def test_cot():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert cot(nan) == nan

    assert cot(oo*I) == -I
    assert cot(-oo*I) == I

    assert cot(0) == cot(0)

    assert cot(1) == cot(1)
    assert cot(-1) == -cot(1)

    assert cot(x) == cot(x)
    assert cot(-x) == -cot(x)

    assert cot(pi*I) == -coth(pi)*I
    assert cot(-pi*I) == coth(pi)*I

    assert cot(2**1024 * E) == cot(2**1024 * E)
    assert cot(-2**1024 * E) == -cot(2**1024 * E)

    assert cot(pi) == cot(pi)
    assert cot(-pi) == -cot(pi)
    assert cot(2*pi) == cot(2*pi)
    assert cot(-2*pi) == -cot(2*pi)
    assert cot(-3*10**73*pi) == -cot(3*10**73*pi)
    assert cot(7*10**103*pi) == cot(7*10**103*pi)

    assert cot(pi/2) == 0
    assert cot(-pi/2) == 0
    assert cot(5*pi/2) == 0
    assert cot(7*pi/2) == 0

    assert cot(pi/3) == 1/sqrt(3)
    assert cot(-2*pi/3) == 1/sqrt(3)

    assert cot(pi/4) == Basic.One()
    assert cot(-pi/4) == -Basic.One()
    assert cot(17*pi/4) == Basic.One()
    assert cot(-3*pi/4) == Basic.One()

    assert cot(pi/6) == sqrt(3)
    assert cot(-pi/6) == -sqrt(3)
    assert cot(7*pi/6) == sqrt(3)
    assert cot(-5*pi/6) == sqrt(3)

    assert cot(pi/105) == cot(pi/105)
    assert cot(-pi/105) == -cot(pi/105)

    assert cot(2 + 3*I) == cot(2 + 3*I)

    assert cot(x*I) == -coth(x)*I

    assert cot(k*pi) == cot(k*pi)
    assert cot(17*k*pi) == cot(17*k*pi)

    assert cot(k*pi*I) == -coth(k*pi)*I

#def test_asin():
#def test_acos():
#def test_atan():
#def test_acot():



def test_attributes():
    x = Symbol('x')
    assert sin(x)[:] == (x,)
    assert sin(x)[0] != sin
    assert sin(x)[0] == x
