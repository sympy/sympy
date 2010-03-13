"""Tests for efficient functions for generating orthogonal polynomials. """

from sympy import Poly, raises, Rational as Q

from sympy.polys.orthopolys import (
    chebyshevt_poly,
    chebyshevu_poly,
    hermite_poly,
    legendre_poly,
    laguerre_poly,
)

from sympy.abc import x

def test_chebyshevt_poly():
    raises(ValueError, "chebyshevt_poly(-1, x)")

    assert chebyshevt_poly(1, x, polys=True) == Poly(x)

    assert chebyshevt_poly(0, x) == 1
    assert chebyshevt_poly(1, x) == x
    assert chebyshevt_poly(2, x) == 2*x**2 - 1
    assert chebyshevt_poly(3, x) == 4*x**3 - 3*x
    assert chebyshevt_poly(4, x) == 8*x**4 - 8*x**2 + 1
    assert chebyshevt_poly(5, x) == 16*x**5 - 20*x**3 + 5*x
    assert chebyshevt_poly(6, x) == 32*x**6 - 48*x**4 + 18*x**2 - 1

def test_chebyshevu_poly():
    raises(ValueError, "chebyshevu_poly(-1, x)")

    assert chebyshevu_poly(1, x, polys=True) == Poly(2*x)

    assert chebyshevu_poly(0, x) == 1
    assert chebyshevu_poly(1, x) == 2*x
    assert chebyshevu_poly(2, x) == 4*x**2 - 1
    assert chebyshevu_poly(3, x) == 8*x**3 - 4*x
    assert chebyshevu_poly(4, x) == 16*x**4 - 12*x**2 + 1
    assert chebyshevu_poly(5, x) == 32*x**5 - 32*x**3 + 6*x
    assert chebyshevu_poly(6, x) == 64*x**6 - 80*x**4 + 24*x**2 - 1

def test_hermite_poly():
    raises(ValueError, "hermite_poly(-1, x)")

    assert hermite_poly(1, x, polys=True) == Poly(2*x)

    assert hermite_poly(0, x) == 1
    assert hermite_poly(1, x) == 2*x
    assert hermite_poly(2, x) == 4*x**2 - 2
    assert hermite_poly(3, x) == 8*x**3 - 12*x
    assert hermite_poly(4, x) == 16*x**4 - 48*x**2 + 12
    assert hermite_poly(5, x) == 32*x**5 - 160*x**3 + 120*x
    assert hermite_poly(6, x) == 64*x**6 - 480*x**4 + 720*x**2 - 120

def test_legendre_poly():
    raises(ValueError, "legendre_poly(-1, x)")

    assert legendre_poly(1, x, polys=True) == Poly(x)

    assert legendre_poly(0, x) == 1
    assert legendre_poly(1, x) == x
    assert legendre_poly(2, x) == Q(3,2)*x**2 - Q(1,2)
    assert legendre_poly(3, x) == Q(5,2)*x**3 - Q(3,2)*x
    assert legendre_poly(4, x) == Q(35,8)*x**4 - Q(30,8)*x**2 + Q(3,8)
    assert legendre_poly(5, x) == Q(63,8)*x**5 - Q(70,8)*x**3 + Q(15,8)*x
    assert legendre_poly(6, x) == Q(231,16)*x**6 - Q(315,16)*x**4 + Q(105,16)*x**2 - Q(5,16)

def test_laguerre_poly():
    raises(ValueError, "laguerre_poly(-1, x)")

    assert laguerre_poly(1, x, polys=True) == Poly(-x + 1)

    assert laguerre_poly(0, x) == 1
    assert laguerre_poly(1, x) == -x + 1
    assert laguerre_poly(2, x) == Q(1,2)*x**2 - Q(4,2)*x + 1
    assert laguerre_poly(3, x) == -Q(1,6)*x**3 + Q(9,6)*x**2 - Q(18,6)*x + 1
    assert laguerre_poly(4, x) == Q(1,24)*x**4 - Q(16,24)*x**3 + Q(72,24)*x**2 - Q(96,24)*x + 1
    assert laguerre_poly(5, x) == -Q(1,120)*x**5 + Q(25,120)*x**4 - Q(200,120)*x**3 + Q(600,120)*x**2 - Q(600,120)*x + 1
    assert laguerre_poly(6, x) == Q(1,720)*x**6 - Q(36,720)*x**5 + Q(450,720)*x**4 - Q(2400,720)*x**3 + Q(5400,720)*x**2 - Q(4320,720)*x + 1

