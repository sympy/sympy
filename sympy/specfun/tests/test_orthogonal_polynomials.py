from sympy import *

x = Symbol('x')
'''
def test_legendre():
    assert legendre(0, x) == 1
    assert legendre(1, x) == x
    assert legendre(2, x) == ((3*x**2-1)/2).expand()
    assert legendre(3, x) == ((5*x**3-3*x)/2).expand()
    assert legendre(10, -1) == 1
    assert legendre(11, -1) == -1
    assert legendre(10, 1) == 1
    assert legendre(11, 1) == 1
    assert legendre(10, 0) != 0
    assert legendre(11, 0) == 0
    """
    for n in range(1, 5):
        for k in range(n):
            z = legendre_zero(n, k)
            assert legendre(n, z) == 0
            assert abs(legendre(n, z.evalf())) < 1e-8
            assert abs(legendre(n+1, z.evalf())) > 1e-8
    """
    assert legendre(3, sqrt(Rational(3,5))) == 0
    assert legendre(3, -sqrt(Rational(3,5))) == 0

def test_chebyshev():
    chebyshev = Chebyshev3()
    assert chebyshev(0, x) == 1
    assert chebyshev(1, x) == x
    assert chebyshev(2, x) == 2*x**2-1
    assert chebyshev(3, x) == 4*x**3-3*x
    """
    for n in range(1, 5):
        for k in range(n):
            z = chebyshev_zero(n, k)
            assert chebyshev(n, z) == 0
            assert abs(chebyshev(n, z.evalf())) < 1e-8
            assert abs(chebyshev(n+1, z.evalf())) > 1e-8
    """
'''
