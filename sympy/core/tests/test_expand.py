from sympy import (log, sqrt, Rational as R, Symbol, I, exp, pi, S, re, im,
    Tuple, sin)

from sympy.simplify.simplify import expand_numer, expand
from sympy.utilities.pytest import raises
from sympy.core.function import (expand_mul, expand_multinomial, expand_log,
    expand_func, expand_trig, expand_complex, expand_power_base,
    expand_power_exp)

from sympy.abc import x, y, z

def test_expand_no_log():
    assert ((1+log(x**4))**2).expand(log=False) == 1 + 2*log(x**4) + log(x**4)**2
    assert ((1+log(x**4))*(1+log(x**3))).expand(log=False) == 1 + log(x**4) + log(x**3) + log(x**4)*log(x**3)

def test_expand_no_multinomial():
    assert ((1+x)*(1+(1+x)**4)).expand(multinomial=False) == 1 + x + (1+x)**4 + x*(1+x)**4

def test_expand_negative_integer_powers():
    expr = (x+y)**(-2)
    assert expr.expand() == 1 / (2*x*y + x**2 + y**2)
    assert expr.expand(multinomial=False) == (x+y)**(-2)
    expr = (x+y)**(-3)
    assert expr.expand() == 1 / (3*x*x*y + 3*x*y*y + x**3 + y**3)
    assert expr.expand(multinomial=False) == (x+y)**(-3)
    expr = (x+y)**(2) * (x+y)**(-4)
    assert expr.expand() == 1 / (2*x*y + x**2 + y**2)
    assert expr.expand(multinomial=False) == (x+y)**(-2)

def test_expand_non_commutative():
    A = x = Symbol('x', commutative=False)
    B = y = Symbol('y', commutative=False)
    a = Symbol('a')
    i = Symbol('i', integer=True)
    assert ((x + y)**2).expand() == x*y + y*x + x**2 + y**2
    assert ((x + y)**3).expand() == (x**2*y + y**2*x + x*y**2 + y*x**2 +
                                     x**3 + y**3 + x*y*x + y*x*y)
    # 3120
    assert ((a*A*B*A**-1)**2).expand() == a**2*A*B**2*A**(-1)
    # Note that (a*A*B*A**-1)**2 is automatically converted to a**2*(A*B*A**-1)**2
    assert ((a*A*B*A**-1)**2).expand(deep=False) == a**2*(A*B*A**-1)**2
    assert ((a*A*B*A**-1)**2).expand(force=True) == a**2*A*B**2*A**(-1)
    assert ((a*A*B)**2).expand() == a**2*A*B*A*B
    assert ((a*A)**2).expand() == a**2*A**2
    assert ((a*A*B)**i).expand() == a**i*(A*B)**i

def test_expand_radicals():
    a = (x + y)**R(1,2)

    assert (a**1).expand() == a
    assert (a**3).expand() == x*a + y*a
    assert (a**5).expand() == x**2*a + 2*x*y*a + y**2*a

    assert (1/a**1).expand() == 1/a
    assert (1/a**3).expand() == 1/(x*a + y*a)
    assert (1/a**5).expand() == 1/(x**2*a + 2*x*y*a + y**2*a)

    a = (x + y)**R(1,3)

    assert (a**1).expand() == a
    assert (a**2).expand() == a**2
    assert (a**4).expand() == x*a + y*a
    assert (a**5).expand() == x*a**2 + y*a**2
    assert (a**7).expand() == x**2*a + 2*x*y*a + y**2*a

def test_expand_modulus():
    assert ((x + y)**11).expand(modulus=11) == x**11 + y**11
    assert ((x + sqrt(2)*y)**11).expand(modulus=11) == x**11 + 10*sqrt(2)*y**11
    assert (x + y/2).expand(modulus=1) == y/2

    raises(ValueError, lambda: ((x + y)**11).expand(modulus=0))
    raises(ValueError, lambda: ((x + y)**11).expand(modulus=x))

def test_issue_2644():
    assert (x*sqrt(x + y)*(1 + sqrt(x + y))).expand() == x**2 + x*y + x*sqrt(x + y)
    assert (x*sqrt(x + y)*(1 + x*sqrt(x + y))).expand() == x**3 + x**2*y + x*sqrt(x + y)

def test_expand_frac():
    assert expand((x + y)*y/x/(x + 1), frac=True) == \
        (x*y + y**2)/(x**2 + x)
    assert expand((x + y)*y/x/(x + 1), numer=True) == \
        (x*y + y**2)/(x*(x + 1))
    assert expand((x + y)*y/x/(x + 1), denom=True) == \
        y*(x + y)/(x**2 + x)
    eq = (x+1)**2/y
    assert expand_numer(eq, multinomial=False) == eq

def test_issue_3022():
    from sympy import cse
    ans = S('''([
        (x0, im(x)),
        (x1, re(x)),
        (x2, atan2(x0, x1)/2),
        (x3, sin(x2)), (x4, cos(x2)),
        (x5, x0**2 + x1**2),
        (x6, atan2(0, x5)/4),
        (x7, cos(x6)),
        (x8, sin(x6)),
        (x9, x4*x7),
        (x10, x4*x8),
        (x11, x3*x8),
        (x12, x3*x7)],
        [sqrt(2)*(x10 + I*x10 + x11 - I*x11 + x12 + I*x12 - x9 + I*x9)/
        (8*pi**(3/2)*x5**(1/4))])''')
    eq = -I*exp(-3*I*pi/4)/(4*pi**(S(3)/2)*sqrt(x))
    r, e = cse((eq).expand(complex=True))
    assert abs((eq - e[0].subs(reversed(r))).subs(x, 1 + 3*I)) < 1e-9

def test_expand_power_base():
    # was test_separate()

    assert expand_power_base((x*y*z)**4) == x**4*y**4*z**4
    assert expand_power_base((x*y*z)**x).is_Pow
    assert expand_power_base((x*y*z)**x, force=True) == x**x*y**x*z**x
    assert expand_power_base((x*(y*z)**2)**3) == x**3*y**6*z**6

    assert expand_power_base((sin((x*y)**2)*y)**z).is_Pow
    assert expand_power_base((sin((x*y)**2)*y)**z, force=True) == sin((x*y)**2)**z*y**z
    assert expand_power_base((sin((x*y)**2)*y)**z, deep=True) == (sin(x**2*y**2)*y)**z

    assert expand_power_base(exp(x)**2) == exp(2*x)
    assert expand_power_base((exp(x)*exp(y))**2) == exp(2*x)*exp(2*y)

    assert expand_power_base((exp((x*y)**z)*exp(y))**2) == exp(2*(x*y)**z)*exp(2*y)
    assert expand_power_base((exp((x*y)**z)*exp(y))**2, deep=True, force=True) == exp(2*x**z*y**z)*exp(2*y)

    assert expand_power_base((exp(x)*exp(y))**z).is_Pow
    assert expand_power_base((exp(x)*exp(y))**z, force=True) == exp(x)**z*exp(y)**z
