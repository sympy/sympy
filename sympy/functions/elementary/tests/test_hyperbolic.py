
from sympy import *

def test_sinh():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert sinh(nan) == nan

    assert sinh(oo) == oo
    assert sinh(-oo) == -oo

    assert sinh(0) == 0

    assert sinh(1) == sinh(1)
    assert sinh(-1) == -sinh(1)

    assert sinh(x) == sinh(x)
    assert sinh(-x) == -sinh(x)

    assert sinh(pi) == sinh(pi)
    assert sinh(-pi) == -sinh(pi)

    assert sinh(2**1024 * E) == sinh(2**1024 * E)
    assert sinh(-2**1024 * E) == -sinh(2**1024 * E)

    assert sinh(pi*I) == 0
    assert sinh(-pi*I) == 0
    assert sinh(2*pi*I) == 0
    assert sinh(-2*pi*I) == 0
    assert sinh(-3*10**73*pi*I) == 0
    assert sinh(7*10**103*pi*I) == 0

    assert sinh(pi*I/2) == I
    assert sinh(-pi*I/2) == -I
    assert sinh(5*pi*I/2) == I
    assert sinh(7*pi*I/2) == -I

    assert sinh(pi*I/3) == Basic.Half()*sqrt(3)*I
    assert sinh(-2*pi*I/3) == -Basic.Half()*sqrt(3)*I

    assert sinh(pi*I/4) == Basic.Half()*sqrt(2)*I
    assert sinh(-pi*I/4) == -Basic.Half()*sqrt(2)*I
    assert sinh(17*pi*I/4) == Basic.Half()*sqrt(2)*I
    assert sinh(-3*pi*I/4) == -Basic.Half()*sqrt(2)*I

    assert sinh(pi*I/6) == Basic.Half()*I
    assert sinh(-pi*I/6) == -Basic.Half()*I
    assert sinh(7*pi*I/6) == -Basic.Half()*I
    assert sinh(-5*pi*I/6) == -Basic.Half()*I

    assert sinh(pi*I/105) == sin(pi/105)*I
    assert sinh(-pi*I/105) == -sin(pi/105)*I

    assert sinh(2 + 3*I) == sinh(2 + 3*I)

    assert sinh(x*I) == sin(x)*I

    assert sinh(k*pi*I) == 0
    assert sinh(17*k*pi*I) == 0

    assert sinh(k*pi*I/2) == sin(k*pi/2)*I

def test_cosh():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert cosh(nan) == nan

    assert cosh(oo) == oo
    assert cosh(-oo) == oo

    assert cosh(0) == 1

    assert cosh(1) == cosh(1)
    assert cosh(-1) == cosh(1)

    assert cosh(x) == cosh(x)
    assert cosh(-x) == cosh(x)

    assert cosh(pi*I) == cos(pi)
    assert cosh(-pi*I) == cos(pi)

    assert cosh(2**1024 * E) == cosh(2**1024 * E)
    assert cosh(-2**1024 * E) == cosh(2**1024 * E)

    assert cosh(pi*I/2) == 0
    assert cosh(-pi*I/2) == 0
    assert cosh(pi*I/2) == 0
    assert cosh(-pi*I/2) == 0
    assert cosh((-3*10**73+1)*pi*I/2) == 0
    assert cosh((7*10**103+1)*pi*I/2) == 0

    assert cosh(pi*I) == -1
    assert cosh(-pi*I) == -1
    assert cosh(5*pi*I) == -1
    assert cosh(8*pi*I) == 1

    assert cosh(pi*I/3) == Basic.Half()
    assert cosh(-2*pi*I/3) == -Basic.Half()

    assert cosh(pi*I/4) == Basic.Half()*sqrt(2)
    assert cosh(-pi*I/4) == Basic.Half()*sqrt(2)
    assert cosh(11*pi*I/4) == -Basic.Half()*sqrt(2)
    assert cosh(-3*pi*I/4) == -Basic.Half()*sqrt(2)

    assert cosh(pi*I/6) == Basic.Half()*sqrt(3)
    assert cosh(-pi*I/6) == Basic.Half()*sqrt(3)
    assert cosh(7*pi*I/6) == -Basic.Half()*sqrt(3)
    assert cosh(-5*pi*I/6) == -Basic.Half()*sqrt(3)

    assert cosh(pi*I/105) == cos(pi/105)
    assert cosh(-pi*I/105) == cos(pi/105)

    assert cosh(2 + 3*I) == cosh(2 + 3*I)

    assert cosh(x*I) == cos(x)

    assert cosh(k*pi*I) == cos(k*pi)
    assert cosh(17*k*pi*I) == cos(17*k*pi)

    assert cosh(k*pi) == cosh(k*pi)

def test_tanh():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert tanh(nan) == nan

    assert tanh(oo) == 1
    assert tanh(-oo) == -1

    assert tanh(0) == 0

    assert tanh(1) == tanh(1)
    assert tanh(-1) == -tanh(1)

    assert tanh(x) == tanh(x)
    assert tanh(-x) == -tanh(x)

    assert tanh(pi) == tanh(pi)
    assert tanh(-pi) == -tanh(pi)

    assert tanh(2**1024 * E) == tanh(2**1024 * E)
    assert tanh(-2**1024 * E) == -tanh(2**1024 * E)

    assert tanh(pi*I) == 0
    assert tanh(-pi*I) == 0
    assert tanh(2*pi*I) == 0
    assert tanh(-2*pi*I) == 0
    assert tanh(-3*10**73*pi*I) == 0
    assert tanh(7*10**103*pi*I) == 0

    assert tanh(pi*I/2) == tanh(pi*I/2)
    assert tanh(-pi*I/2) == -tanh(pi*I/2)
    assert tanh(5*pi*I/2) == tanh(5*pi*I/2)
    assert tanh(7*pi*I/2) == tanh(7*pi*I/2)

    assert tanh(pi*I/3) == sqrt(3)*I
    assert tanh(-2*pi*I/3) == sqrt(3)*I

    assert tanh(pi*I/4) == I
    assert tanh(-pi*I/4) == -I
    assert tanh(17*pi*I/4) == I
    assert tanh(-3*pi*I/4) == I

    assert tanh(pi*I/6) == I/sqrt(3)
    assert tanh(-pi*I/6) == -I/sqrt(3)
    assert tanh(7*pi*I/6) == I/sqrt(3)
    assert tanh(-5*pi*I/6) == I/sqrt(3)

    assert tanh(pi*I/105) == tan(pi/105)*I
    assert tanh(-pi*I/105) == -tan(pi/105)*I

    assert tanh(2 + 3*I) == tanh(2 + 3*I)

    assert tanh(x*I) == tan(x)*I

    assert tanh(k*pi*I) == 0
    assert tanh(17*k*pi*I) == 0

    assert tanh(k*pi*I/2) == tan(k*pi/2)*I

def test_coth():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert coth(nan) == nan

    assert coth(oo) == 1
    assert coth(-oo) == -1

    assert coth(0) == coth(0)

    assert coth(1) == coth(1)
    assert coth(-1) == -coth(1)

    assert coth(x) == coth(x)
    assert coth(-x) == -coth(x)

    assert coth(pi*I) == -cot(pi)*I
    assert coth(-pi*I) == cot(pi)*I

    assert coth(2**1024 * E) == coth(2**1024 * E)
    assert coth(-2**1024 * E) == -coth(2**1024 * E)

    assert coth(pi*I) == -cot(pi)*I
    assert coth(-pi*I) == cot(pi)*I
    assert coth(2*pi*I) == -cot(2*pi)*I
    assert coth(-2*pi*I) == cot(2*pi)*I
    assert coth(-3*10**73*pi*I) == cot(3*10**73*pi)*I
    assert coth(7*10**103*pi*I) == -cot(7*10**103*pi)*I

    assert coth(pi*I/2) == 0
    assert coth(-pi*I/2) == 0
    assert coth(5*pi*I/2) == 0
    assert coth(7*pi*I/2) == 0

    assert coth(pi*I/3) == -I/sqrt(3)
    assert coth(-2*pi*I/3) == -I/sqrt(3)

    assert coth(pi*I/4) == -I
    assert coth(-pi*I/4) == I
    assert coth(17*pi*I/4) == -I
    assert coth(-3*pi*I/4) == -I

    assert coth(pi*I/6) == -sqrt(3)*I
    assert coth(-pi*I/6) == sqrt(3)*I
    assert coth(7*pi*I/6) == -sqrt(3)*I
    assert coth(-5*pi*I/6) == -sqrt(3)*I

    assert coth(pi*I/105) == -cot(pi/105)*I
    assert coth(-pi*I/105) == cot(pi/105)*I

    assert coth(2 + 3*I) == coth(2 + 3*I)

    assert coth(x*I) == -cot(x)*I

    assert coth(k*pi*I) == -cot(k*pi)*I
    assert coth(17*k*pi*I) == -cot(17*k*pi)*I

    assert coth(k*pi*I) == -cot(k*pi)*I

#def test_asinh():
#def test_acosh():
#def test_atanh():
#def test_acoth():
