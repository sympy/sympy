from sympy import symbols, Symbol, nan, oo, I, sinh, sin, acot, pi, atan, \
        acos, Rational, sqrt, asin, acot, cot, coth, E, S, tan, tanh, cos, \
        cosh, atan2, exp

def test_sin():
    x, y = symbols('xy')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert sin(nan) == nan

    assert sin(oo*I) == oo*I
    assert sin(-oo*I) == -oo*I

    assert sin(0) == 0

    assert sin(1) == sin(1)
    assert sin(-1) == -sin(1)

    assert sin(x) == sin(x)
    assert sin(-x) == -sin(x)

    assert sin(asin(x)) == x
    assert sin(atan(x)) == x / sqrt(1 + x**2)
    assert sin(acos(x)) == sqrt(1 - x**2)
    assert sin(acot(x)) == 1 / (sqrt(1 + 1 / x**2) * x)

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

    assert sin(pi/3) == S.Half*sqrt(3)
    assert sin(-2*pi/3) == -S.Half*sqrt(3)

    assert sin(pi/4) == S.Half*sqrt(2)
    assert sin(-pi/4) == -S.Half*sqrt(2)
    assert sin(17*pi/4) == S.Half*sqrt(2)
    assert sin(-3*pi/4) == -S.Half*sqrt(2)

    assert sin(pi/6) == S.Half
    assert sin(-pi/6) == -S.Half
    assert sin(7*pi/6) == -S.Half
    assert sin(-5*pi/6) == -S.Half

    assert sin(1*pi/5) == sqrt((5 - sqrt(5)) / 8)
    assert sin(2*pi/5) == sqrt((5 + sqrt(5)) / 8)
    assert sin(3*pi/5) == sin(2*pi/5)
    assert sin(4*pi/5) == sin(1*pi/5)
    assert sin(6*pi/5) == -sin(1*pi/5)
    assert sin(8*pi/5) == -sin(2*pi/5)

    assert sin(-1273*pi/5) == -sin(2*pi/5)

    assert sin(pi/105) == sin(pi/105)
    assert sin(-pi/105) == -sin(pi/105)

    assert sin(104*pi/105) == sin(pi/105)
    assert sin(106*pi/105) == -sin(pi/105)

    assert sin(-104*pi/105) == -sin(pi/105)
    assert sin(-106*pi/105) == sin(pi/105)

    assert sin(2 + 3*I) == sin(2 + 3*I)

    assert sin(x*I) == sinh(x)*I

    assert sin(k*pi) == 0
    assert sin(17*k*pi) == 0

    assert sin(k*pi*I) == sinh(k*pi)*I

    assert sin(r).is_real == True

    assert sin(exp(10)-1) == sin(-1+exp(10))

def test_cos():
    x, y = symbols('xy')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert cos(nan) == nan

    assert cos(oo*I) == oo
    assert cos(-oo*I) == oo

    assert cos(0) == 1

    assert cos(1) == cos(1)
    assert cos(-1) == cos(1)

    assert cos(x) == cos(x)
    assert cos(-x) == cos(x)

    assert cos(acos(x)) == x
    assert cos(atan(x)) == 1 / sqrt(1 + x**2)
    assert cos(asin(x)) == sqrt(1 - x**2)
    assert cos(acot(x)) == 1 / sqrt(1 + 1 / x**2)

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
    assert cos(2*pi)==1
    assert cos(5*pi) == -1
    assert cos(8*pi) == 1

    assert cos(pi/3) == S.Half
    assert cos(-2*pi/3) == -S.Half

    assert cos(pi/4) == S.Half*sqrt(2)
    assert cos(-pi/4) == S.Half*sqrt(2)
    assert cos(11*pi/4) == -S.Half*sqrt(2)
    assert cos(-3*pi/4) == -S.Half*sqrt(2)

    assert cos(pi/6) == S.Half*sqrt(3)
    assert cos(-pi/6) == S.Half*sqrt(3)
    assert cos(7*pi/6) == -S.Half*sqrt(3)
    assert cos(-5*pi/6) == -S.Half*sqrt(3)

    assert cos(1*pi/5) == (sqrt(5) + 1)/4
    assert cos(2*pi/5) == (sqrt(5) - 1)/4
    assert cos(3*pi/5) == -cos(2*pi/5)
    assert cos(4*pi/5) == -cos(1*pi/5)
    assert cos(6*pi/5) == -cos(1*pi/5)
    assert cos(8*pi/5) == cos(2*pi/5)

    assert cos(-1273*pi/5) == -cos(2*pi/5)

    assert cos(pi/105) == cos(pi/105)
    assert cos(-pi/105) == cos(pi/105)

    assert cos(104*pi/105) == -cos(pi/105)
    assert cos(106*pi/105) == -cos(pi/105)

    assert cos(-104*pi/105) == -cos(pi/105)
    assert cos(-106*pi/105) == -cos(pi/105)

    assert cos(2 + 3*I) == cos(2 + 3*I)

    assert cos(x*I) == cosh(x)

    assert cos(k*pi) == cos(k*pi)
    assert cos(17*k*pi) == cos(17*k*pi)

    assert cos(k*pi*I) == cosh(k*pi)

    assert cos(r).is_real == True

    assert cos(exp(10)-1) == cos(-1+exp(10))

def test_tan():
    x, y = symbols('xy')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert tan(nan) == nan

    assert tan(oo*I) == I
    assert tan(-oo*I) == -I

    assert tan(0) == 0

    assert tan(1) == tan(1)
    assert tan(-1) == -tan(1)

    assert tan(x) == tan(x)
    assert tan(-x) == -tan(x)

    assert tan(atan(x)) == x
    assert tan(asin(x)) == x / sqrt(1 - x**2)
    assert tan(acos(x)) == sqrt(1 - x**2) / x
    assert tan(acot(x)) == 1 / x

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

    assert tan(pi/4) == S.One
    assert tan(-pi/4) == -S.One
    assert tan(17*pi/4) == S.One
    assert tan(-3*pi/4) == S.One

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

    assert tan(r).is_real == True

def test_cot():
    x, y = symbols('xy')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    assert cot(nan) == nan

    assert cot(oo*I) == -I
    assert cot(-oo*I) == I

    assert cot(0) == cot(0)

    assert cot(1) == cot(1)
    assert cot(-1) == -cot(1)

    assert cot(x) == cot(x)
    assert cot(-x) == -cot(x)

    assert cot(acot(x)) == x
    assert cot(atan(x)) == 1 / x
    assert cot(asin(x)) == sqrt(1 - x**2) / x
    assert cot(acos(x)) == x / sqrt(1 - x**2)

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

    assert cot(pi/4) == S.One
    assert cot(-pi/4) == -S.One
    assert cot(17*pi/4) == S.One
    assert cot(-3*pi/4) == S.One

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

    assert cot(r).is_real == True

# TODO write me
def test_asin():
    x = Symbol('x')

    assert asin(0)  == 0
    assert asin(Rational(1,2)) == pi/6
    assert asin(1)  == pi/2
    assert asin(sqrt(3)/2) == pi/3

    assert asin(x).diff(x) ==  1/sqrt(1-x**2)

    assert asin(0.2).is_real == True
    assert asin(-2).is_real == False

# TODO write me
def test_acos():
    x = Symbol('x')

    r = Symbol('r', real=True)

    assert acos(0)  == pi/2
    assert acos(Rational(1,2)) == pi/3
    assert acos(1)  == 0
    assert acos(sqrt(2)/2) == pi/4
    assert acos(x).diff(x) == -1/sqrt(1-x**2)

    assert acos(0.2).is_real == True
    assert acos(-2).is_real == False


# TODO write me
def test_atan():
    x = Symbol('x')

    r = Symbol('r', real=True)

    assert atan(0)  == 0
    assert atan(1)  == pi/4
    assert atan(sqrt(3)) == pi/3
    assert atan(oo) == pi/2
    assert atan(x).diff(x) ==  1/(1+x**2)

    assert atan(r).is_real == True

def test_atan2():
    assert atan2(0, 0) == S.NaN
    assert atan2(0, 1) == 0
    assert atan2(1, 0) == pi/2
    assert atan2(1, -1) == 3*pi/4
    assert atan2(-1, 1) == -pi/4
    assert atan2(0, -1) == pi

# TODO write me
def test_acot():
    x = Symbol('x')

    r = Symbol('r', real=True)

    assert acot(oo) == 0
    assert acot(1)  == pi/4
    assert acot(0)  == pi/2
    assert acot(x).diff(x) == -1/(1+x**2)

    assert acot(r).is_real == True



def test_attributes():
    x = Symbol('x')
    assert sin(x).args[:] == (x,)
    assert sin(x).args[0] != sin
    assert sin(x).args[0] == x

def test_sincos_rewrite():
    x = Symbol("x")
    y = Symbol("y")
    assert sin(pi/2-x) == cos(x)
    assert sin(pi-x) == sin(x)
    assert cos(pi/2-x) == sin(x)
    assert cos(pi-x) == -cos(x)

    assert sin(-x-y) == -sin(x+y)
    assert sin(-x-1) == -sin(x+1)
    assert cos(-x-y) == cos(x+y)
    assert cos(-x-1) == cos(x+1)
    assert sin(x-y) == sin(x-y)
    assert sin(y-x) == sin(y-x)
    assert cos(y-x) == cos(y-x)
    assert cos(x-y) == cos(x-y)

