from sympy.sandbox.core import Symbol, Rational, sqrt

a = Symbol("a")
b = Symbol("b")
c = Symbol("c")

def test_Symbol():
    assert str(a)=="a"
    assert str(b)=="b"
    e=a*b
    assert e==a*b
    assert a*b*b==a*b**2
    assert a*b*b+c==c+a*b**2
    assert a*b*b-c==-c+a*b**2

def test_arit():
    p = Rational(5)
    e=a*b
    assert e == a*b
    e=a*b+b*a
    assert e == 2*a*b,`e,2*a*b`
    e=a*b+b*a+a*b+p*b*a
    assert e == 8*a*b
    e=a*b+b*a+a*b+p*b*a+a
    assert e == a+8*a*b
    e=a+a
    assert e == 2*a
    e=a+b+a
    assert e == b+2*a
    e=a+b*b+a+b*b
    assert e == 2*a+2*b**2
    e=a+Rational(2)+b*b+a+b*b+p
    assert e == 7+2*a+2*b**2
    e=(a+b*b+a+b*b)*p
    assert e == 5*(2*a+2*b**2)
    e=(a*b*c+c*b*a+b*a*c)*p
    assert e == 15*a*b*c
    e=(a*b*c+c*b*a+b*a*c)*p-Rational(15)*a*b*c
    assert e == Rational(0)
    e = Rational(50)*(a-a)
    assert e == Rational(0)
    e=b*a-b-a*b+b
    assert e == Rational(0)
    e=a*b+c**p
    assert e == a*b+c**5
    e=a/b
    assert e == a*b**(-1)
    e=a*2*2
    assert e == 4*a
    e=2+a*2/2
    assert e == 2+a
    e=2-a-2
    assert e == -a
    e=2*a*2
    assert e == 4*a
    e=2/a/2
    assert e == a**(-1)
    e=2**a**2
    assert e == 2**(a**2)
    e = -(1+a)
    assert e == -1 -a
    e = Rational(1,2)*(1+a)
    assert e == Rational(1,2) + a/2

def testdiv():
    e=a/b
    assert e == a*b**(-1)
    e=a/b+c/2
    assert e == a*b**(-1)+Rational(1)/2*c
    e=(1-b)/(b-1)
    assert e == (1+-b)*((-1)+b)**(-1)

def testpow():
    n1 = Rational(1)
    n2 = Rational(2)
    n5 = Rational(5)
    e=a*a
    assert e == a**2
    e=a*a*a
    assert e == a**3
    e=a*a*a*a**Rational(6)
    assert e == a**9
    e=a*a*a*a**Rational(6)-a**Rational(9)
    assert e == Rational(0)
    e=a**(b+c)*a**(-b)
    assert e == a**c
    e=a**(b+c)*a*a**(-b)*a**(-c)/a
    assert e == Rational(1)
    e=a**(b-b)
    assert e == Rational(1)
    e=(a-a)**b
    #assert e == Rational(0),`e` # this is true only when real part of b is positive
    e=(a+Rational(1)-a)**b
    assert e == Rational(1)

    e=(a+b+c)**n2
    assert e == (a+b+c)**2
    #assert e.expand() == 2*b*c+2*a*c+2*a*b+a**2+c**2+b**2

    e=(a+b)**n2
    assert e == (a+b)**2
    #assert e.expand() == 2*a*b+a**2+b**2

    e=(a+b)**(n1/n2)
    assert e == (a+b)**(Rational(1)/2)
    #assert e.expand() == (a+b)**(Rational(1)/2)

    n=n5**(n1/n2)
    assert n == Rational(5)**(Rational(1)/2)
    e=n*a*b-n*b*a
    assert e == Rational(0)
    e=n*a*b+n*b*a
    assert e == 2*a*b*5**(Rational(1)/2)
    #assert e.diff(a) == 2*b*5**(Rational(1)/2)
    #assert e.diff(a) == 2*b*5**(Rational(1)/2)
    e=a/b**2
    assert e == a*b**(-2)

    assert sqrt(2*(1+sqrt(2))) == (2*(1+2**(Rational(1,2))))**(Rational(1,2))

    x = Symbol('x')
    y = Symbol('y')

    #assert ((x*y)**3).expand() == y**3 * x**3
    #assert ((x*y)**-3).expand() == y**-3 * x**-3

    #assert (x**5*(3*x)**(3)).expand() == 27 * x**8
    #assert (x**5*(-3*x)**(3)).expand() == -27 * x**8
    #assert (x**5*(3*x)**(-3)).expand() == Rational(1,27) * x**2
    #assert (x**5*(-3*x)**(-3)).expand() == -Rational(1,27) * x**2

    #assert ((-x)**2)**Rational(1,3) == ((-x)**Rational(1,3))**2
    #assert (-x)**Rational(2,3) == x**Rational(2,3),`(-x)**Rational(2,3), x**Rational(2,3)`
    #assert (-x)**Rational(5,7) == -x**Rational(5,7)
    #assert 4**Rational(1, 4) == 2**Rational(1, 2)

def test_expand():
    p = Rational(5)
    e = (a+b)*c
    assert e == c*(a+b)
    #assert (e.expand()-a*c-b*c) == Rational(0)
    e=(a+b)*(a+b)
    assert e == (a+b)**2
    #assert e.expand() == 2*a*b+a**2+b**2
    e=(a+b)*(a+b)**Rational(2)
    assert e == (a+b)**3
    #assert e.expand() == 3*b*a**2+3*a*b**2+a**3+b**3
    #assert e.expand() == 3*b*a**2+3*a*b**2+a**3+b**3
    e=(a+b)*(a+c)*(b+c)
    assert e == (a+c)*(a+b)*(b+c)
    #assert e.expand() == 2*a*b*c+b*a**2+c*a**2+b*c**2+a*c**2+c*b**2+a*b**2
    e=(a+Rational(1))**p
    assert e == (1+a)**5
    #assert e.expand() == 1+5*a+10*a**2+10*a**3+5*a**4+a**5
    e=(a+b+c)*(a+c+p)
    assert e == (5+a+c)*(a+b+c)
    #assert e.expand() == 5*a+5*b+5*c+2*a*c+b*c+a*b+a**2+c**2
    #x=Symbol("x")
    #s=exp(x*x)-1
    #e=s.series(x,5)/x**2
    #assert e.expand() ==  1+x**2/2+O(x**3)

    # Check that this isn't too slow
    #x = Symbol('x')
    #W = 1
    #for i in range(1, 21):
    #    W = W * (x-i)
    #W = W.expand()
    #assert W.has(-1672280820*x**15)


def test_power_expand():
    """Test for Pow.expand()"""
    a = Symbol('a')
    b = Symbol('b')
    p = (a+b)**2
    #assert p.expand() == a**2 + b**2 + 2*a*b

    p = (1+2*(1+a))**2
    #assert p.expand() == 9 + 4*(a**2) + 12*a

def test_powerbug():
    x=Symbol("x")
    assert x**1 != (-x)**1
    assert x**2 == (-x)**2
    assert x**3 != (-x)**3
    assert x**4 == (-x)**4
    assert x**5 != (-x)**5
    assert x**6 == (-x)**6

    assert x**128 == (-x)**128
    assert x**129 != (-x)**129

    assert (2*x)**2 == (-2*x)**2

if __name__ == "__main__":
    # make output from recursion errors more pleasant
    #import sys
    #sys.setrecursionlimit(30)
    test_Symbol()
    test_arit()
    testdiv()
    testpow()
    test_expand()
    test_power_expand()
    test_powerbug()
    print "passed"
