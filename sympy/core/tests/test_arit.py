from sympy import Symbol, sin, cos, exp, O, sqrt, Rational

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

def _test_arit():
    p = Rational(5)
    e=a*b
    assert e == a*b
    e=a*b+b*a
    assert e == 2*a*b
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

def _testpow():
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
    # XXX Fails, returns 0**b instead of 0
    e=(a-a)**b
    assert e == Rational(0)
    e=(a+Rational(1)-a)**b
    assert e == Rational(1)

    e=(a+b+c)**n2
    assert e == (a+b+c)**2
    assert e.expand() == 2*b*c+2*a*c+2*a*b+a**2+c**2+b**2

    e=(a+b)**n2
    assert e == (a+b)**2
    assert e.expand() == 2*a*b+a**2+b**2

    e=(a+b)**(n1/n2)
    assert e == (a+b)**(Rational(1)/2)
    assert e.expand() == (a+b)**(Rational(1)/2)

    n=n5**(n1/n2)
    assert n == Rational(5)**(Rational(1)/2)
    e=n*a*b-n*b*a
    assert e == Rational(0)
    e=n*a*b+n*b*a
    assert e == 2*a*b*5**(Rational(1)/2)
    assert e.diff(a) == 2*b*5**(Rational(1)/2)
    assert e.diff(a) == 2*b*5**(Rational(1)/2)
    e=a/b**2
    assert e == a*b**(-2)

    assert sqrt(2*(1+sqrt(2))) == (2*(1+2**(Rational(1,2))))**(Rational(1,2))

    x = Symbol('x')
    y = Symbol('y')

    assert ((x*y)**3).expand() == y**3 * x**3
    assert ((x*y)**-3).expand() == y**-3 * x**-3

    assert (x**5*(3*x)**(3)).expand() == 27 * x**8
    assert (x**5*(-3*x)**(3)).expand() == -27 * x**8
    assert (x**5*(3*x)**(-3)).expand() == Rational(1,27) * x**2
    assert (x**5*(-3*x)**(-3)).expand() == -Rational(1,27) * x**2

    n = Symbol('k', even=False)
    k = Symbol('k', even=True)

    assert (-1)**x == (-1)**x
    assert (-1)**n == (-1)**n
    assert (-2)**k == 2**k
    assert (-1)**k == 1

    assert ((-x)**2)**Rational(1,3) == ((-x)**Rational(1,3))**2
    assert (-x)**Rational(2,3) == x**Rational(2,3)
    assert (-x)**Rational(5,7) == -x**Rational(5,7)
    assert 4**Rational(1, 4) == 2**Rational(1, 2)

def test_expand():
    p = Rational(5)
    e = (a+b)*c
    assert e == c*(a+b)
    assert (e.expand()-a*c-b*c) == Rational(0)
    e=(a+b)*(a+b)
    assert e == (a+b)**2
    assert e.expand() == 2*a*b+a**2+b**2
    e=(a+b)*(a+b)**Rational(2)
    assert e == (a+b)**3
    assert e.expand() == 3*b*a**2+3*a*b**2+a**3+b**3
    assert e.expand() == 3*b*a**2+3*a*b**2+a**3+b**3
    e=(a+b)*(a+c)*(b+c)
    assert e == (a+c)*(a+b)*(b+c)
    assert e.expand() == 2*a*b*c+b*a**2+c*a**2+b*c**2+a*c**2+c*b**2+a*b**2
    e=(a+Rational(1))**p
    assert e == (1+a)**5
    assert e.expand() == 1+5*a+10*a**2+10*a**3+5*a**4+a**5
    e=(a+b+c)*(a+c+p)
    assert e == (5+a+c)*(a+b+c)
    assert e.expand() == 5*a+5*b+5*c+2*a*c+b*c+a*b+a**2+c**2
    x=Symbol("x")
    s=exp(x*x)-1
    e=s.series(x,5)/x**2
    #assert e == (x**2+x**4/2)/x**2
    assert e.expand() ==  1+x**2/2+O(x**3)

def test_power_expand():
    """Test for Pow.expand()"""
    a = Symbol('a')
    b = Symbol('b')
    p = (a+b)**2
    assert p.expand() == a**2 + b**2 + 2*a*b

    p = (1+2*(1+a))**2
    assert p.expand() == 9 + 4*(a**2) + 12*a

def test_ncmul():
    A = Symbol("A", commutative=False)
    B = Symbol("B", commutative=False)
    C = Symbol("C", commutative=False)
    assert A*B != B*A
    assert A*B*C != C*B*A
    assert A*b*B*3*C == 3*b*A*B*C
    assert A*b*B*3*C != 3*b*B*A*C
    assert A*b*B*3*C == 3*A*B*C*b

    assert A+B == B+A
    assert (A+B)*C != C*(A+B)

    assert C*(A+B)*C != C*C*(A+B)

    assert (C*(A+B)).expand() == C*A+C*B
    assert (C*(A+B)).expand() != A*C+B*C

    assert A*A == A**2
    assert (A+B)*(A+B) == (A+B)**2
    assert ((A+B)**2).expand() == A**2 + A*B + B*A +B**2

    assert A**-1  * A == 1
    assert A/A == 1
    assert A/(A**2) == 1/A

    assert A/(1+A) == A/(1+A)

def test_ncpow():
    x = Symbol('x', commutative=False)
    y = Symbol('y', commutative=False)

    assert (x**2)*(y**2) != (y**2)*(x**2)
    assert (x**-2)*y != y*(x**2)

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

def test_Add_Mul_is_integer():
    x = Symbol('x')

    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True)

    assert (2*k).is_integer == True
    assert (-k).is_integer == True
    assert (k/3).is_integer == False
    assert (x*k*n).is_integer == None

    assert (k+n).is_integer == True
    assert (k+x).is_integer == None
    assert (k+n*x).is_integer == None
    assert (k+n/3).is_integer == False

def _test_Add_Mul_is_bounded():
    x = Symbol('x')

    # XXX Fails, and the assertion isn't correct. The assertion below can
    #     only be made when x is real.
    assert sin(x).is_bounded == True
    assert (x*sin(x)).is_bounded == None
    assert (1024*sin(x)).is_bounded == True
    assert (sin(x)*exp(x)).is_bounded == False
    assert (sin(x)*cos(x)).is_bounded == True
    assert (x*sin(x)*exp(x)).is_bounded == False

    assert (sin(x)-67).is_bounded == True
    assert (sin(x)+exp(x)).is_bounded == False

def test_Mul_is_even_odd():
    x = Symbol('x', integer=True)

    k = Symbol('k', odd=True)
    n = Symbol('n', odd=True)
    m = Symbol('m', even=True)

    assert (2*x).is_even == True
    assert (2*x).is_odd == False

    assert (3*x).is_even == None
    assert (3*x).is_odd == None

    assert (k/3).is_even == None
    assert (k/3).is_odd == None

    assert (2*n).is_even == True
    assert (2*n).is_odd == False

    assert (2*m).is_even == True
    assert (2*m).is_odd == False

    assert (-n).is_even == False
    assert (-n).is_odd == True

    assert (k*n).is_even == False
    assert (k*n).is_odd == True

    assert (k*m).is_even == True
    assert (k*m).is_odd == False

    assert (k*n*m).is_even == True
    assert (k*n*m).is_odd == False

    assert (k*m*x).is_even == True
    assert (k*m*x).is_odd == False

def test_Add_is_even_odd():
    x = Symbol('x', integer=True)

    k = Symbol('k', odd=True)
    n = Symbol('n', even=True)

    assert (2+k).is_even == False
    assert (2+k).is_odd == True

    assert (7-k).is_even == True
    assert (7-k).is_odd == False

    assert (11-n).is_even == False
    assert (11-n).is_odd == True

    assert (-8+n).is_even == True
    assert (-8+n).is_odd == False

    assert (n+k).is_even == False
    assert (n+k).is_odd == True

    assert (n-k).is_even == False
    assert (n-k).is_odd == True

    assert (n+2*k).is_even == True
    assert (n+2*k).is_odd == False

    assert (k+n+x).is_odd == None
    assert (k+n-x).is_even == None

    assert (2*k+n*x).is_odd == None
    assert (2*k+n*x).is_even == None

def _test_Mul_is_negative_positive():
    x = Symbol('x', real=True)
    y = Symbol('y', real=False)

    k = Symbol('k', negative=True)
    n = Symbol('n', positive=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)

    assert k.is_negative == True
    assert (-k).is_negative == False
    assert (2*k).is_negative == True

    assert n.is_negative == False
    assert (-n).is_negative == True
    # XXX Fails, but the assertion is correct
    assert (2*n).is_negative == False

    assert (n*k).is_negative == True
    assert (2*n*k).is_negative == True
    assert (-n*k).is_negative == False
    assert (n*k*y).is_negative == None

    assert u.is_negative == False
    assert (-u).is_negative == None
    assert (2*u).is_negative == False

    assert v.is_negative == None
    assert (-v).is_negative == False
    assert (2*v).is_negative == None

    assert (u*v).is_negative == None

    assert (k*u).is_negative == None
    assert (k*v).is_negative == False

    assert (n*u).is_negative == False
    assert (n*v).is_negative == None

    assert (v*k*u).is_negative == False
    assert (v*n*u).is_negative == None

    assert (-v*k*u).is_negative == None
    assert (-v*n*u).is_negative == False

    assert (17*v*k*u).is_negative == False
    assert (17*v*n*u).is_negative == None

    assert (k*v*n*u).is_negative == False

    assert (x*k).is_negative == None
    assert (u*v*n*x*k).is_negative == None

    assert k.is_positive == False
    assert (-k).is_positive == True
    assert (2*k).is_positive == False

    assert n.is_positive == True
    assert (-n).is_positive == False
    assert (2*n).is_positive == True

    assert (n*k).is_positive == False
    assert (2*n*k).is_positive == False
    assert (-n*k).is_positive == True
    assert (-n*k*y).is_positive == None

    assert u.is_positive == None
    assert (-u).is_positive == False
    assert (2*u).is_positive == None

    assert v.is_positive == False
    assert (-v).is_positive == None
    assert (2*v).is_positive == False

    assert (u*v).is_positive == False

    assert (k*u).is_positive == False
    assert (k*v).is_positive == None

    assert (n*u).is_positive == None
    assert (n*v).is_positive == False

    assert (v*k*u).is_positive == None
    assert (v*n*u).is_positive == False

    assert (-v*k*u).is_positive == False
    assert (-v*n*u).is_positive == None

    assert (17*v*k*u).is_positive == None
    assert (17*v*n*u).is_positive == False

    assert (k*v*n*u).is_positive == None

    assert (x*k).is_positive == None
    assert (u*v*n*x*k).is_positive == None

def _test_Mul_is_nonpositive_nonnegative():
    x = Symbol('x', real=True)

    k = Symbol('k', negative=True)
    n = Symbol('n', positive=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)

    # XXX Fails, but the assertion is correct
    assert k.is_nonpositive == True
    assert (-k).is_nonpositive == False
    assert (2*k).is_nonpositive == True

    assert n.is_nonpositive == False
    assert (-n).is_nonpositive == True
    assert (2*n).is_nonpositive == False

    assert (n*k).is_nonpositive == True
    assert (2*n*k).is_nonpositive == True
    assert (-n*k).is_nonpositive == False

    assert u.is_nonpositive == None
    assert (-u).is_nonpositive == True
    assert (2*u).is_nonpositive == None

    assert v.is_nonpositive == True
    assert (-v).is_nonpositive == None
    assert (2*v).is_nonpositive == True

    assert (u*v).is_nonpositive == True

    assert (k*u).is_nonpositive == True
    assert (k*v).is_nonpositive == None

    assert (n*u).is_nonpositive == None
    assert (n*v).is_nonpositive == True

    assert (v*k*u).is_nonpositive == None
    assert (v*n*u).is_nonpositive == True

    assert (-v*k*u).is_nonpositive == True
    assert (-v*n*u).is_nonpositive == None

    assert (17*v*k*u).is_nonpositive == None
    assert (17*v*n*u).is_nonpositive == True

    assert (k*v*n*u).is_nonpositive == None

    assert (x*k).is_nonpositive == None
    assert (u*v*n*x*k).is_nonpositive == None

    assert k.is_nonnegative == False
    assert (-k).is_nonnegative == True
    assert (2*k).is_nonnegative == False

    assert n.is_nonnegative == True
    assert (-n).is_nonnegative == False
    assert (2*n).is_nonnegative == True

    assert (n*k).is_nonnegative == False
    assert (2*n*k).is_nonnegative == False
    assert (-n*k).is_nonnegative == True

    assert u.is_nonnegative == True
    assert (-u).is_nonnegative == None
    assert (2*u).is_nonnegative == True

    assert v.is_nonnegative == None
    assert (-v).is_nonnegative == True
    assert (2*v).is_nonnegative == None

    assert (u*v).is_nonnegative == None

    assert (k*u).is_nonnegative == None
    assert (k*v).is_nonnegative == True

    assert (n*u).is_nonnegative == True
    assert (n*v).is_nonnegative == None

    assert (v*k*u).is_nonnegative == True
    assert (v*n*u).is_nonnegative == None

    assert (-v*k*u).is_nonnegative == None
    assert (-v*n*u).is_nonnegative == True

    assert (17*v*k*u).is_nonnegative == True
    assert (17*v*n*u).is_nonnegative == None

    assert (k*v*n*u).is_nonnegative == True

    assert (x*k).is_nonnegative == None
    assert (u*v*n*x*k).is_nonnegative == None

def test_Add_is_even_odd():
    x = Symbol('x', integer=True)

    k = Symbol('k', odd=True)
    n = Symbol('n', odd=True)
    m = Symbol('m', even=True)

    assert (k+7).is_even == True
    assert (k+7).is_odd == False

    assert (-k+7).is_even == True
    assert (-k+7).is_odd == False

    assert (k-12).is_even == False
    assert (k-12).is_odd == True

    assert (-k-12).is_even == False
    assert (-k-12).is_odd == True

    assert (k+n).is_even == True
    assert (k+n).is_odd == False

    assert (k+m).is_even == False
    assert (k+m).is_odd == True

    assert (k+n+m).is_even == True
    assert (k+n+m).is_odd == False

    assert (k+n+x+m).is_even == None
    assert (k+n+x+m).is_odd == None

def _test_Add_is_negative_positive():
    x = Symbol('x', real=True)

    k = Symbol('k', negative=True)
    n = Symbol('n', positive=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)

    assert (k-2).is_negative == True
    assert (k+17).is_negative == None
    assert (-k-5).is_negative == None
    assert (-k+123).is_negative == False

    assert (k-n).is_negative == True
    assert (k+n).is_negative == None
    assert (-k-n).is_negative == None
    # XXX Fails, but the assertion is correct
    assert (-k+n).is_negative == False

    assert (k-n-2).is_negative == True
    assert (k+n+17).is_negative == None
    assert (-k-n-5).is_negative == None
    assert (-k+n+123).is_negative == False

    assert (-2*k+123*n+17).is_negative == False

    assert (k+u).is_negative == None
    assert (k+v).is_negative == None
    assert (n+u).is_negative == False
    assert (n+v).is_negative == None

    assert (u-v).is_negative == False
    assert (u+v).is_negative == None
    assert (-u-v).is_negative == None
    assert (-u+v).is_negative == None

    assert (u-v+n+2).is_negative == False
    assert (u+v+n+2).is_negative == None
    assert (-u-v+n+2).is_negative == None
    assert (-u+v+n+2).is_negative == None

    assert (k+x).is_negative == None
    assert (k+x-n).is_negative == None

    assert (k-2).is_positive == False
    assert (k+17).is_positive == None
    assert (-k-5).is_positive == None
    assert (-k+123).is_positive == True

    assert (k-n).is_positive == False
    assert (k+n).is_positive == None
    assert (-k-n).is_positive == None
    assert (-k+n).is_positive == True

    assert (k-n-2).is_positive == False
    assert (k+n+17).is_positive == None
    assert (-k-n-5).is_positive == None
    assert (-k+n+123).is_positive == True

    assert (-2*k+123*n+17).is_positive == True

    assert (k+u).is_positive == None
    assert (k+v).is_positive == False
    assert (n+u).is_positive == None
    assert (n+v).is_positive == None

    assert (u-v).is_positive == None
    assert (u+v).is_positive == None
    assert (-u-v).is_positive == None
    assert (-u+v).is_positive == False

    assert (u-v-n-2).is_positive == None
    assert (u+v-n-2).is_positive == None
    assert (-u-v-n-2).is_positive == None
    assert (-u+v-n-2).is_positive == False

    assert (n+x).is_positive == None
    assert (n+x-k).is_positive == None

def _test_Add_is_nonpositive_nonnegative():
    x = Symbol('x', real=True)

    k = Symbol('k', negative=True)
    n = Symbol('n', positive=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)

    assert (u-2).is_nonpositive == None
    assert (u+17).is_nonpositive == False
    assert (-u-5).is_nonpositive == True
    assert (-u+123).is_nonpositive == None

    assert (u-v).is_nonpositive == None
    assert (u+v).is_nonpositive == None
    assert (-u-v).is_nonpositive == None
    assert (-u+v).is_nonpositive == True

    assert (u-v-2).is_nonpositive == None
    assert (u+v+17).is_nonpositive == None
    assert (-u-v-5).is_nonpositive == None
    assert (-u+v-123).is_nonpositive == True

    assert (-2*u+123*v-17).is_nonpositive == True

    assert (k+u).is_nonpositive == None
    # XXX Fails, but the assertion is correct
    assert (k+v).is_nonpositive == True
    assert (n+u).is_nonpositive == None
    assert (n+v).is_nonpositive == None

    assert (k-n).is_nonpositive == True
    assert (k+n).is_nonpositive == None
    assert (-k-n).is_nonpositive == None
    assert (-k+n).is_nonpositive == False

    assert (k-n+u+2).is_nonpositive == None
    assert (k+n+u+2).is_nonpositive == None
    assert (-k-n+u+2).is_nonpositive == None
    assert (-k+n+u+2).is_nonpositive == None

    assert (u+x).is_nonpositive == None
    assert (v-x-n).is_nonpositive == None

    assert (u-2).is_nonnegative == None
    assert (u+17).is_nonnegative == True
    assert (-u-5).is_nonnegative == None
    assert (-u+123).is_nonnegative == None

    assert (u-v).is_nonnegative == True
    assert (u+v).is_nonnegative == None
    assert (-u-v).is_nonnegative == None
    assert (-u+v).is_nonnegative == None

    assert (u-v+2).is_nonnegative == True
    assert (u+v+17).is_nonnegative == None
    assert (-u-v-5).is_nonnegative == None
    assert (-u+v-123).is_nonnegative == None

    assert (2*u-123*v+17).is_nonnegative == True

    assert (k+u).is_nonnegative == None
    assert (k+v).is_nonnegative == None
    assert (n+u).is_nonnegative == True
    assert (n+v).is_nonnegative == None

    assert (k-n).is_nonnegative == False
    assert (k+n).is_nonnegative == None
    assert (-k-n).is_nonnegative == None
    assert (-k+n).is_nonnegative == True

    assert (k-n-u-2).is_nonnegative == None
    assert (k+n-u-2).is_nonnegative == None
    assert (-k-n-u-2).is_nonnegative == None
    assert (-k+n-u-2).is_nonnegative == None

    assert (u-x).is_nonnegative == None
    assert (v+x+n).is_nonnegative == None

def _test_Pow_is_integer():
    x = Symbol('x')

    k = Symbol('k', integer=True)
    n = Symbol('n', nni=True)
    m = Symbol('m', pi=True)

    assert (k**2).is_integer == True
    assert (k**(-2)).is_integer == False

    assert (2**k).is_integer == None
    assert (2**(-k)).is_integer == None

    assert (2**n).is_integer == True
    # XXX Fails, returns False (n could be 0, and hence this value
    #     would still be an integer, so the assertion is correct)
    assert (2**(-n)).is_integer == None

    assert (2**m).is_integer == True
    assert (2**(-m)).is_integer == False

    assert (x**2).is_integer == None
    assert (2**x).is_integer == None

    assert (k**n).is_integer == True
    assert (k**(-n)).is_integer == None

    assert (k**x).is_integer == None
    assert (x**k).is_integer == None

    assert (k**(n*m)).is_integer == True
    assert (k**(-n*m)).is_integer == None

def _test_Pow_is_bounded():
    x = Symbol('x')

    # XXX is failing for some reason
    assert (x**2).is_bounded == None

    assert (sin(x)**2).is_bounded == True
    assert (sin(x)**x).is_bounded == None
    assert (sin(x)**exp(x)).is_bounded == False

    assert (1/sin(x)).is_bounded == False
    assert (1/exp(x)).is_bounded == False

def _test_Pow_is_even_odd():
    x = Symbol('x')

    k = Symbol('k', even=True)
    n = Symbol('n', odd=True)
    m = Symbol('m', nni=True)

    assert (k**2).is_even == True
    assert (n**2).is_even == False
    assert (2**k).is_even == None

    assert (k**m).is_even == True
    assert (n**m).is_even == False

    assert (k**x).is_even == None
    # XXX Fails
    assert (n**x).is_even == False

    assert (k**2).is_odd == False
    assert (n**2).is_odd == True
    assert (3**k).is_odd == None

    assert (k**m).is_odd == False
    assert (n**m).is_odd == True

    assert (k**x).is_odd == False
    assert (n**x).is_odd == None

def test_Pow_is_negative_positive():
    x = Symbol('x', real=True)

    k = Symbol('k', pi=True)
    n = Symbol('n', even=True)
    m = Symbol('m', odd=True)

    assert (2**x).is_positive == True
    assert ((-2)**x).is_positive == None
    assert ((-2)**n).is_positive == True
    assert ((-2)**m).is_positive == False

    assert (k**2).is_positive == True
    assert (k**(-2)).is_positive == True

    assert (k**x).is_positive == True
    assert ((-k)**x).is_positive == None
    assert ((-k)**n).is_positive == True
    assert ((-k)**m).is_positive == False

    assert (2**x).is_negative == False
    assert ((-2)**x).is_negative == None
    assert ((-2)**n).is_negative == False
    assert ((-2)**m).is_negative == True

    assert (k**2).is_negative == False
    assert (k**(-2)).is_negative == False

    assert (k**x).is_negative == False
    assert ((-k)**x).is_negative == None
    assert ((-k)**n).is_negative == False
    assert ((-k)**m).is_negative == True

def test_Pow_is_nonpositive_nonnegative():
    x = Symbol('x', real=True)

    k = Symbol('k', nni=True)
    n = Symbol('n', even=True)
    m = Symbol('m', odd=True)

    assert (2**x).is_nonnegative == True
    assert ((-2)**x).is_nonnegative == None
    assert ((-2)**n).is_nonnegative == True
    assert ((-2)**m).is_nonnegative == False

    assert (k**2).is_nonnegative == True
    assert (k**(-2)).is_nonnegative == True

    assert (k**x).is_nonnegative == True
    assert ((-k)**x).is_nonnegative == None
    assert ((-k)**n).is_nonnegative == True
    assert ((-k)**m).is_nonnegative == False

    assert (2**x).is_nonpositive == False
    assert ((-2)**x).is_nonpositive == None
    assert ((-2)**n).is_nonpositive == False
    assert ((-2)**m).is_nonpositive == True

    assert (k**2).is_nonpositive == None
    assert (k**(-2)).is_nonpositive == None

    assert (k**x).is_nonpositive == None
    assert ((-k)**x).is_nonpositive == None
    assert ((-k)**n).is_nonpositive == False
    assert ((-k)**m).is_nonpositive == True
