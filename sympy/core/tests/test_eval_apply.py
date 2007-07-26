
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

    x, y = symbols('xy', real=True)

    assert log(x) == log(x)
    assert log(x*y) == log(x) + log(y)

    assert log(x**2) == 2*log(x)
    assert log(x**y) == log(x**y)

    k = Symbol('k', positive=True)

    assert log(-x) == log(-x)
    assert log(-k) == log(-k)

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

def test_apply_evalf():
    value = (log(3)/log(2)-1).evalf()

    assert value.epsilon_eq(Real("0.58496250072115618145373"))
