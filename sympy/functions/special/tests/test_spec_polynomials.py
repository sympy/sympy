from sympy import (legendre, Symbol, Dummy, diff, Derivative, Rational, roots, sympify, S, sqrt,
                   cos, pi, binomial, Sum, hermite, chebyshevu, chebyshevt, chebyshevt_root,
                   chebyshevu_root, assoc_legendre, laguerre, assoc_laguerre, laguerre_poly, sqrt)

from sympy.utilities.pytest import raises

x = Symbol('x')

def test_legendre():
    raises(ValueError, lambda: legendre(-1, x))
    assert legendre(0, x) == 1
    assert legendre(1, x) == x
    assert legendre(2, x) == ((3*x**2-1)/2).expand()
    assert legendre(3, x) == ((5*x**3-3*x)/2).expand()
    assert legendre(4, x) == ((35*x**4-30*x**2+3)/8).expand()
    assert legendre(5, x) == ((63*x**5-70*x**3+15*x)/8).expand()
    assert legendre(6, x) == ((231*x**6-315*x**4+105*x**2-5)/16).expand()

    assert legendre(10, -1) == 1
    assert legendre(11, -1) == -1
    assert legendre(10, 1) == 1
    assert legendre(11, 1) == 1
    assert legendre(10, 0) != 0
    assert legendre(11, 0) == 0

    assert roots(legendre(4,x), x) == {
         sqrt(Rational(3, 7) - Rational(2, 35)*sqrt(30)): 1,
        -sqrt(Rational(3, 7) - Rational(2, 35)*sqrt(30)): 1,
         sqrt(Rational(3, 7) + Rational(2, 35)*sqrt(30)): 1,
        -sqrt(Rational(3, 7) + Rational(2, 35)*sqrt(30)): 1,
    }

    n = Symbol("n")

    X = legendre(n,x)
    assert isinstance(X, legendre)

    assert legendre(-n,x) == legendre(n-1, x)
    assert legendre(n,-x) == (-1)**n*legendre(n, x)
    assert diff(legendre(n,x), x) == n*(x*legendre(n, x) - legendre(n - 1, x))/(x**2 - 1)
    assert diff(legendre(n,x), n) == Derivative(legendre(n, x), n)

def test_assoc_legendre():
    Plm=assoc_legendre
    Q=sqrt(1-x**2)

    assert Plm(0, 0, x) ==  1
    assert Plm(1, 0, x) ==  x
    assert Plm(1, 1, x) == -Q
    assert Plm(2, 0, x) ==  (3*x**2-1)/2
    assert Plm(2, 1, x) == -3*x*Q
    assert Plm(2, 2, x) ==  3*Q**2
    assert Plm(3, 0, x) ==  (5*x**3-3*x)/2
    assert Plm(3, 1, x).expand() ==  (( 3*(1-5*x**2)/2 ).expand() * Q).expand()
    assert Plm(3, 2, x) ==  15*x * Q**2
    assert Plm(3, 3, x) == -15 * Q**3

    # negative m
    assert Plm(1,-1, x) == -Plm(1, 1, x)/2
    assert Plm(2,-2, x) ==  Plm(2, 2, x)/24
    assert Plm(2,-1, x) == -Plm(2, 1, x)/6
    assert Plm(3,-3, x) == -Plm(3, 3, x)/720
    assert Plm(3,-2, x) ==  Plm(3, 2, x)/120
    assert Plm(3,-1, x) == -Plm(3, 1, x)/12

    n = Symbol("n")
    m = Symbol("m")

    X = Plm(n,m, x)
    assert isinstance(X, assoc_legendre)

    assert Plm(n,0, x) == legendre(n, x)

    raises(ValueError, lambda: Plm(-1, 0, x))
    raises(ValueError, lambda: Plm(0, 1, x))

def test_chebyshev():

    assert chebyshevt(0, x) == 1
    assert chebyshevt(1, x) == x
    assert chebyshevt(2, x) == 2*x**2-1
    assert chebyshevt(3, x) == 4*x**3-3*x

    for n in range(1, 4):
        for k in range(n):
            z = chebyshevt_root(n, k)
            assert chebyshevt(n, z) == 0
        raises(ValueError, lambda: chebyshevt_root(n, n))

    for n in range(1, 4):
        for k in range(n):
            z = chebyshevu_root(n, k)
            assert chebyshevu(n, z) == 0
        raises(ValueError, lambda: chebyshevu_root(n, n))

    n = Symbol("n")
    X = chebyshevt(n, x)
    assert isinstance(X, chebyshevt)
    assert chebyshevt(n, -x) == (-1)**n*chebyshevt(n, x)
    assert chebyshevt(-n, x) == chebyshevt(n, x)

    assert chebyshevt(n, 0) == cos(pi*n/2)
    assert chebyshevt(n, 1) == 1

    assert diff(chebyshevt(n, x), x) == n*chebyshevu(n - 1, x)

    X = chebyshevu(n, x)
    assert isinstance(X, chebyshevu)

    assert chebyshevu(n, -x) == (-1)**n*chebyshevu(n, x)
    assert chebyshevu(-n, x) == -chebyshevu(n - 2, x)

    assert chebyshevu(n, 0) == cos(pi*n/2)
    assert chebyshevu(n, 1) == n + 1

    assert diff(chebyshevu(n, x), x) == (-x*chebyshevu(n, x) + (n + 1)*chebyshevt(n + 1, x))/(x**2 - 1)

def test_hermite():
    assert hermite(0, x) == 1
    assert hermite(1, x) == 2*x
    assert hermite(2, x) == 4*x**2 - 2
    assert hermite(3, x) == 8*x**3 - 12*x
    assert hermite(4, x) == 16*x**4 - 48*x**2 + 12
    assert hermite(6, x) == 64*x**6 - 480*x**4 + 720*x**2 - 120

    n = Symbol("n")
    assert hermite(n, x) == hermite(n, x)
    assert hermite(n, -x) == (-1)**n*hermite(n, x)
    assert hermite(-n, x) == hermite(-n, x)

    assert diff(hermite(n, x), x) == 2*n*hermite(n - 1, x)
    assert diff(hermite(n, x), n) ==  Derivative(hermite(n, x), n)

def test_laguerre():
    alpha = Symbol("alpha")

    # generalized Laguerre polynomials:
    assert assoc_laguerre(0, alpha, x) == 1
    assert assoc_laguerre(1, alpha, x) == -x + alpha + 1
    assert assoc_laguerre(2, alpha, x).expand() == (x**2/2 - (alpha+2)*x + (alpha+2)*(alpha+1)/2).expand()
    assert assoc_laguerre(3, alpha, x).expand() == (-x**3/6 + (alpha+3)*x**2/2 - (alpha+2)*(alpha+3)*x/2 + (alpha+1)*(alpha+2)*(alpha+3)/6).expand()

    # Laguerre polynomials:
    assert assoc_laguerre(0, 0, x) == 1
    assert assoc_laguerre(1, 0, x) == 1 - x
    assert assoc_laguerre(2, 0, x).expand() == 1 - 2*x + x**2/2
    assert assoc_laguerre(3, 0, x).expand() == 1 - 3*x + 3*x**2/2 - x**3/6

    # Test the lowest 10 polynomials with laguerre_poly, to make sure that it
    # works:
    for i in range(10):
        assert assoc_laguerre(i, 0, x).expand() == laguerre_poly(i, x)

    n = Symbol("n")
    X = laguerre(n, x)
    assert isinstance(X, laguerre)

    assert laguerre(n, 0) == 1

    assert diff(laguerre(n, x), x) == -assoc_laguerre(n - 1, 1, x)

    m = Symbol("m")
    X = assoc_laguerre(n, m, x)
    assert isinstance(X, assoc_laguerre)

    assert assoc_laguerre(n, 0, x) == laguerre(n, x)
    assert assoc_laguerre(n, alpha, 0) == binomial(alpha + n, alpha)

    assert diff(assoc_laguerre(n, alpha, x), x) == -assoc_laguerre(n - 1, alpha + 1, x)
    #k = Dummy("k")
    #assert diff(assoc_laguerre(n, alpha, x), alpha) == Sum(assoc_laguerre(k, alpha, x)/(-alpha + n), (k, 0, n - 1))
