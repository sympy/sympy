from sympy.core import *

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

    assert feq(Basic.One() / Real(-0.5), -Integer(2))

def test_expand():
    x = Symbol('x')
    assert (2**(-1-x)).expand() == Rational(1,2)*2**(-x)

def test_issue153():
    #test that is runs:
    a = Basic.sqrt(2*(1+Basic.sqrt(2)))

def test_issue350():
    #test if powers are simplified correctly
    a = Symbol('a')
    assert ((a**Rational(1,3))**Rational(2)) == a**Rational(2,3)
    assert ((a**Rational(3))**Rational(2,5)) != a**Rational(6,5)

    a = Symbol('a', real = True)
    assert (a**Rational(3))**Rational(2,5) == a**Rational(6,5)
    
    #assert Number(5)**Rational(2,3)==Number(25)**Rational(1,3)

    
test_issue350()
