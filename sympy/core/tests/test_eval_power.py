from sympy.core import Rational, Symbol, Basic, S, Float, Integer
from sympy.functions.elementary.miscellaneous import sqrt

def test_rational():
    a = Rational(1, 5)

    assert a**Rational(1, 2) == a**Rational(1, 2)
    assert 2 * a**Rational(1, 2) == 2 * a**Rational(1, 2)

    assert a**Rational(3, 2) == a * a**Rational(1, 2)
    assert 2 * a**Rational(3, 2) == 2*a * a**Rational(1, 2)

    assert a**Rational(17, 3) == a**5 * a**Rational(2, 3)
    assert 2 * a**Rational(17, 3) == 2*a**5 * a**Rational(2, 3)

def test_large_rational():
    e = (Rational(123712**12-1,7)+Rational(1,7))**Rational(1,3)
    assert e == 234232585392159195136 * (Rational(1,7)**Rational(1,3))

def test_negative_real():
    def feq(a,b):
        return abs(a - b) < 1E-10

    assert feq(S.One / Float(-0.5), -Integer(2))

def test_expand():
    x = Symbol('x')
    assert (2**(-1-x)).expand() == Rational(1,2)*2**(-x)

def test_issue153():
    #test that is runs:
    a = sqrt(2*(1+sqrt(2)))

def test_issue350():
    #test if powers are simplified correctly
    #see also issue 896
    x = Symbol('x')
    assert ((x**Rational(1,3))**Rational(2)) == x**Rational(2,3)
    assert ((x**Rational(3))**Rational(2,5)) == (x**Rational(3))**Rational(2,5)

    a = Symbol('a', real=True)
    b = Symbol('b', real=True)
    assert (a**2)**b == (abs(a)**b)**2
    assert sqrt(1/a) != 1/sqrt(a) # e.g. for a = -1
    assert (a**3)**Rational(1,3) != a
    assert (x**a)**b != x**(a*b) # e.g. x = -1, a=1/2, b=2
    assert (x**.5)**b == x**(.5*b)
    assert (x**.5)**.5 == x**.25
    assert (x**2.5)**.5 != x**1.25 # e.g. for x = 5*I

    k = Symbol('k',integer=True)
    m = Symbol('m',integer=True)
    assert (x**k)**m == x**(k*m)
    #assert Number(5)**Rational(2,3)==Number(25)**Rational(1,3)

    assert (x**.5)**2 == x**1.0
    assert (x**2)**k == (x**k)**2 == x**(2*k)

    a = Symbol('a', positive=True)
    assert (a**3)**Rational(2,5) == a**Rational(6,5)
    assert (a**2)**b == (a**b)**2
    assert (a**Rational(2, 3))**x == (a**x)**Rational(2, 3) == (a**(2*x/3))

def test_issue767():
    assert --sqrt(sqrt(5)-1)==sqrt(sqrt(5)-1)

def test_negative_one():
    x = Symbol('x', complex=True)
    y = Symbol('y', complex=True)
    assert 1/x**y == x**(-y)

def test_issue1263():
    neg = Symbol('neg', negative=True)
    nonneg = Symbol('nonneg', negative=False)
    any = Symbol('any')
    num, den = sqrt(1/neg).as_numer_denom()
    assert num == sqrt(-1)
    assert den == sqrt(-neg)
    num, den = sqrt(1/nonneg).as_numer_denom()
    assert num == 1
    assert den == sqrt(nonneg)
    num, den = sqrt(1/any).as_numer_denom()
    assert num == sqrt(1/any)
    assert den == 1

    def eqn(num, den, pow):
        return (num/den)**pow
    npos=1
    nneg=-1
    dpos=2-sqrt(3)
    dneg=1-sqrt(3)
    I = S.ImaginaryUnit
    assert dpos > 0 and dneg < 0 and npos > 0 and nneg < 0
    # pos or neg integer
    eq=eqn(npos, dpos, 2);assert eq.is_Pow and eq.as_numer_denom() == (1, dpos**2)
    eq=eqn(npos, dneg, 2);assert eq.is_Pow and eq.as_numer_denom() == (1, dneg**2)
    eq=eqn(nneg, dpos, 2);assert eq.is_Pow and eq.as_numer_denom() == (1, dpos**2)
    eq=eqn(nneg, dneg, 2);assert eq.is_Pow and eq.as_numer_denom() == (1, dneg**2)
    eq=eqn(npos, dpos, -2);assert eq.is_Pow and eq.as_numer_denom() == (dpos**2, 1)
    eq=eqn(npos, dneg, -2);assert eq.is_Pow and eq.as_numer_denom() == (dneg**2, 1)
    eq=eqn(nneg, dpos, -2);assert eq.is_Pow and eq.as_numer_denom() == (dpos**2, 1)
    eq=eqn(nneg, dneg, -2);assert eq.is_Pow and eq.as_numer_denom() == (dneg**2, 1)
    # pos or neg rational
    pow = S.Half
    eq=eqn(npos, dpos, pow);assert eq.is_Pow and eq.as_numer_denom() == (npos**pow, dpos**pow)
    eq=eqn(npos, dneg, pow);assert eq.is_Pow and eq.as_numer_denom() == ((-npos)**pow, (-dneg)**pow)
    eq=eqn(nneg, dpos, pow);assert not eq.is_Pow or eq.as_numer_denom() == (nneg**pow, dpos**pow)
    eq=eqn(nneg, dneg, pow);assert eq.is_Pow and eq.as_numer_denom() == ((-nneg)**pow, (-dneg)**pow)
    eq=eqn(npos, dpos, -pow);assert eq.is_Pow and eq.as_numer_denom() == (dpos**pow, npos**pow)
    eq=eqn(npos, dneg, -pow);assert eq.is_Pow and eq.as_numer_denom() == ((-dneg)**pow, (-npos)**pow)
    eq=eqn(nneg, dpos, -pow);assert not eq.is_Pow or eq.as_numer_denom() == (dpos**pow, nneg**pow)
    eq=eqn(nneg, dneg, -pow);assert eq.is_Pow and eq.as_numer_denom() == ((-dneg)**pow, (-nneg)**pow)
    # unknown exponent
    pow = 2*any
    eq=eqn(npos, dpos, pow)
    assert eq.is_Pow and eq.as_numer_denom() == (npos**pow, dpos**pow)
    eq=eqn(npos, dneg, pow)
    assert eq.is_Pow and eq.as_numer_denom() == (eq, 1)
    eq=eqn(nneg, dpos, pow)
    assert eq.is_Pow and eq.as_numer_denom() == (nneg**pow, dpos**pow)
    eq=eqn(nneg, dneg, pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-nneg)**pow, (-dneg)**pow)
    eq=eqn(npos, dpos, -pow)
    assert eq.as_numer_denom() == (dpos**pow, npos**pow)
    eq=eqn(npos, dneg, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == (1, eq.base**pow)
    eq=eqn(nneg, dpos, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == (dpos**pow, nneg**pow)
    eq=eqn(nneg, dneg, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-dneg)**pow, (-nneg)**pow)

    x = Symbol('x')
    assert ((1/(1 + x/3))**(-S.One)).as_numer_denom() == (3 + x, 3)
    np = Symbol('np',positive=False)
    assert (((1 + x/np)**-2)**(-S.One)).as_numer_denom() == ((np + x)**2, np**2)

def test_Pow_signs():
    """Cf. issues 1496 and 2151"""
    x = Symbol('x')
    y = Symbol('y')
    n = Symbol('n', even=True)
    assert (3-y)**2 != (y-3)**2
    assert (3-y)**n != (y-3)**n
    assert (-3+y-x)**2 != (3-y+x)**2
    assert (y-3)**3 != -(3-y)**3

def test_power_with_noncommutative_mul_as_base():
    x = Symbol('x', commutative=False)
    y = Symbol('y', commutative=False)
    assert not (x*y)**3 == x**3*y**3
    assert (2*x*y)**3 == 8*(x*y)**3

def test_zero():
    x = Symbol('x')
    y = Symbol('y')
    assert 0**x != 0
    assert 0**(2*x) == 0**x
    assert (0**(2 - x)).as_base_exp() == (0, 2 - x)
    assert 0**(x - 2) != S.Infinity**(2 - x)
    assert 0**(2*x*y) == 0**(x*y)
    assert 0**(-2*x*y) == S.Infinity**(x*y)

def test_pow_as_base_exp():
    x = Symbol('x')
    assert (S.Infinity**(2 - x)).as_base_exp() == (S.Infinity, 2 - x)
    assert (S.Infinity**(x - 2)).as_base_exp() == (S.Infinity, x - 2)
