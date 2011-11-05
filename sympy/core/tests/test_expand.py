from sympy import log, sqrt, Rational as R, Symbol

from sympy.simplify.simplify import expand_numer, expand
from sympy.utilities.pytest import raises
from sympy.abc import x, y

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

def test_expand_non_commutative_multinomial():
    x = Symbol('x', commutative=False)
    y = Symbol('x', commutative=False)
    assert ((x+y)**2).expand() == x*y + y*x + x**2 + y**2
    assert ((x+y)**3).expand() == x**2*y + y**2*x + x*y**2 + y*x**2 + x**3 + y**3 + x*y*x + y*x*y

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

    raises(ValueError, "((x + y)**11).expand(modulus=0)")
    raises(ValueError, "((x + y)**11).expand(modulus=x)")

def test_issue_2644():
    assert (x*sqrt(x + y)*(1 + sqrt(x + y))).expand() == x**2 + x*y + x*sqrt(x + y)

def test_expand_frac():
    assert expand((x + y)*y/x/(x + 1), frac=True) == \
        (x*y + y**2)/(x**2 + x)
    assert expand((x + y)*y/x/(x + 1), numer=True) == \
        (x*y + y**2)/(x*(x + 1))
    assert expand((x + y)*y/x/(x + 1), denom=True) == \
        y*(x + y)/(x**2 + x)
    eq = (x+1)**2/y
    assert expand_numer(eq, multinomial=False) == eq
