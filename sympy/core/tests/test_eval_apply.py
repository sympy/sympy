
from sympy import *

def test_log_apply_eval():

    assert log(nan) == nan

    assert log(oo) == oo
    assert log(-oo) == oo

    assert log(0) == -oo

    assert log(1) == 0
    assert log(-1) == I*pi

    assert log(E) == 1
    assert log(-E) == 1 + I*pi

    assert log(pi) == log(pi)
    assert log(-pi) == log(pi) + I*pi

    assert log(17) == log(17)
    assert log(-17) == log(17) + I*pi

    assert log(I) == I*pi/2
    assert log(-I) == -I*pi/2

    assert log(17*I) == I*pi/2 + log(17)
    assert log(-17*I) == -I*pi/2 + log(17)

    assert log(oo*I) == oo
    assert log(-oo*I) == oo

    assert exp(-log(3))**(-1) == 3

    x, y = symbols('xy')

    assert log(x) == log(x)
    assert log(x*y) == log(x*y)

    assert log(x**2) == 2*log(x)
    assert log(x**y) == log(x**y)

    assert log(exp(x)) == log(exp(x))

    x, y = symbols('xy', real=True)

    assert log(x) == log(x)
    assert log(x*y) == log(x) + log(y)

    assert log(x**2) == 2*log(x)
    assert log(x**y) == log(x**y)

    assert log(exp(x)) == x
    assert log(-exp(x)) == x + I*pi

    k = Symbol('k', positive=True)

    assert log(-x) == log(-x)
    assert log(-k) == log(-k)

def test_log_apply_evalf():
    value = (log(3)/log(2)-1).evalf()

    assert value.epsilon_eq(Real("0.58496250072115618145373"))

def test_floor_apply_eval():

    x = Symbol('x')
    y = Symbol('y', real=True)
    k, n = symbols('kn', integer=True)

    assert floor(nan) == nan

    assert floor(oo) == oo
    assert floor(-oo) == -oo

    assert floor(0) == 0

    assert floor(1) == 1
    assert floor(-1) == -1

    assert floor(E) == 2
    assert floor(-E) == -3

    assert floor(pi) == 3
    assert floor(-pi) == -4

    assert floor(Rational(1, 2)) == 0
    assert floor(-Rational(1, 2)) == -1

    assert floor(Rational(7, 3)) == 2
    assert floor(-Rational(7, 3)) == -3

    assert floor(Real(17.0)) == 17
    assert floor(-Real(17.0)) == -17

    assert floor(Real(7.69)) == 7
    assert floor(-Real(7.69)) == -8

    assert floor(I) == I
    assert floor(-I) == -I

    assert floor(oo*I) == oo*I
    assert floor(-oo*I) == -oo*I

    assert floor(2*I) == 2*I
    assert floor(-2*I) == -2*I

    assert floor(I/2) == 0
    assert floor(-I/2) == -I

    assert floor(E + 17) == 19
    assert floor(pi + 2) == 5

    assert floor(E + pi) == floor(E + pi)
    assert floor(I + pi) == floor(I + pi)

    assert floor(floor(pi)) == 3
    assert floor(floor(y)) == floor(y)
    assert floor(floor(x)) == floor(floor(x))

    assert floor(x) == floor(x)
    assert floor(2*x) == floor(2*x)
    assert floor(k*x) == floor(k*x)

    assert floor(k) == k
    assert floor(2*k) == 2*k
    assert floor(k*n) == k*n

    assert floor(k/2) == floor(k/2)

    assert floor(x + y) == floor(x + y)
    assert floor(x + k) == floor(x) + k
    assert floor(k + n) == k + n

    assert floor(x*I) == floor(x*I)
    assert floor(k*I) == k*I

    assert floor(Rational(23, 10) - E*I) == 2 - 3*I

    assert floor(sin(1)) == 0
    assert floor(sin(-1)) == -1

def test_ceiling_apply_eval():

    x = Symbol('x')
    y = Symbol('y', real=True)
    k, n = symbols('kn', integer=True)

    assert ceiling(nan) == nan

    assert ceiling(oo) == oo
    assert ceiling(-oo) == -oo

    assert ceiling(0) == 0

    assert ceiling(1) == 1
    assert ceiling(-1) == -1

    assert ceiling(E) == 3
    assert ceiling(-E) == -2

    assert ceiling(pi) == 4
    assert ceiling(-pi) == -3

    assert ceiling(Rational(1, 2)) == 1
    assert ceiling(-Rational(1, 2)) == 0

    assert ceiling(Rational(7, 3)) == 3
    assert ceiling(-Rational(7, 3)) == -2

    assert ceiling(Real(17.0)) == 17
    assert ceiling(-Real(17.0)) == -17

    assert ceiling(Real(7.69)) == 8
    assert ceiling(-Real(7.69)) == -7

    assert ceiling(I) == I
    assert ceiling(-I) == -I

    assert ceiling(oo*I) == oo*I
    assert ceiling(-oo*I) == -oo*I

    assert ceiling(2*I) == 2*I
    assert ceiling(-2*I) == -2*I

    assert ceiling(I/2) == I
    assert ceiling(-I/2) == 0

    assert ceiling(E + 17) == 20
    assert ceiling(pi + 2) == 6

    assert ceiling(E + pi) == ceiling(E + pi)
    assert ceiling(I + pi) == ceiling(I + pi)

    assert ceiling(ceiling(pi)) == 4
    assert ceiling(ceiling(y)) == ceiling(y)
    assert ceiling(ceiling(x)) == ceiling(ceiling(x))

    assert ceiling(x) == ceiling(x)
    assert ceiling(2*x) == ceiling(2*x)
    assert ceiling(k*x) == ceiling(k*x)

    assert ceiling(k) == k
    assert ceiling(2*k) == 2*k
    assert ceiling(k*n) == k*n

    assert ceiling(k/2) == ceiling(k/2)

    assert ceiling(x + y) == ceiling(x + y)
    assert ceiling(x + k) == ceiling(x) + k
    assert ceiling(k + n) == k + n

    assert ceiling(x*I) == ceiling(x*I)
    assert ceiling(k*I) == k*I

    assert ceiling(Rational(23, 10) - E*I) == 3 - 2*I

    assert ceiling(sin(1)) == 1
    assert ceiling(sin(-1)) == 0

def test_rf_apply_eval():

    x, y = symbols('xy')

    assert rf(nan, y) == nan

    assert rf(x, y) == rf(x, y)

    assert rf(oo, 0) == 1
    assert rf(-oo, 0) == 1

    assert rf(oo, 6) == oo
    assert rf(-oo, 7) == -oo

    assert rf(oo, -6) == oo
    assert rf(-oo, -7) == oo

    assert rf(x, 0) == 1
    assert rf(x, 1) == x
    assert rf(x, 2) == x*(x+1)
    assert rf(x, 3) == x*(x+1)*(x+2)

    assert rf(x, -1) == 1/(x-1)
    assert rf(x, -2) == 1/((x-1)*(x-2))
    assert rf(x, -3) == 1/((x-1)*(x-2)*(x-3))

    assert rf(1, 100) == factorial(100)

def test_ff_apply_eval():

    x, y = symbols('xy')

    assert ff(nan, y) == nan

    assert ff(x, y) == ff(x, y)

    assert ff(oo, 0) == 1
    assert ff(-oo, 0) == 1

    assert ff(oo, 6) == oo
    assert ff(-oo, 7) == -oo

    assert ff(oo, -6) == oo
    assert ff(-oo, -7) == oo

    assert ff(x, 0) == 1
    assert ff(x, 1) == x
    assert ff(x, 2) == x*(x-1)
    assert ff(x, 3) == x*(x-1)*(x-2)

    assert ff(x, -1) == 1/(x+1)
    assert ff(x, -2) == 1/((x+1)*(x+2))
    assert ff(x, -3) == 1/((x+1)*(x+2)*(x+3))

    assert ff(100, 100) == factorial(100)

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
