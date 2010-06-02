"""Tests for functions for generating interesting polynomials. """

from sympy import Poly, ZZ
from sympy.utilities.pytest import raises

from sympy.polys.specialpolys import (
    swinnerton_dyer_poly,
    cyclotomic_poly,
    symmetric_poly,
    random_poly,
    fateman_poly_F_1,
    dmp_fateman_poly_F_1,
    fateman_poly_F_2,
    dmp_fateman_poly_F_2,
    fateman_poly_F_3,
    dmp_fateman_poly_F_3,
)

from sympy.utilities import all, any

from sympy.abc import x, y, z

def test_swinnerton_dyer_poly():
    raises(ValueError, "swinnerton_dyer_poly(0, x)")

    assert swinnerton_dyer_poly(1, x, polys=True) == Poly(x**2 - 2)

    assert swinnerton_dyer_poly(1, x) == x**2 - 2
    assert swinnerton_dyer_poly(2, x) == x**4 - 10*x**2 + 1
    assert swinnerton_dyer_poly(3, x) == x**8 - 40*x**6 + 352*x**4 - 960*x**2 + 576

def test_cyclotomic_poly():
    raises(ValueError, "cyclotomic_poly(0, x)")

    assert cyclotomic_poly(1, x, polys=True) == Poly(x - 1)

    assert cyclotomic_poly(1, x) == x - 1
    assert cyclotomic_poly(2, x) == x + 1
    assert cyclotomic_poly(3, x) == x**2 + x + 1
    assert cyclotomic_poly(4, x) == x**2 + 1
    assert cyclotomic_poly(5, x) == x**4 + x**3 + x**2 + x + 1
    assert cyclotomic_poly(6, x) == x**2 - x + 1

def test_symmetric_poly():
    raises(ValueError, "symmetric_poly(-1, x, y, z)")
    raises(ValueError, "symmetric_poly(5, x, y, z)")

    assert symmetric_poly(1, x, y, z, polys=True) == Poly(x + y + z)
    assert symmetric_poly(1, (x, y, z), polys=True) == Poly(x + y + z)

    assert symmetric_poly(0, x, y, z) == 1
    assert symmetric_poly(1, x, y, z) == x + y + z
    assert symmetric_poly(2, x, y, z) == x*y + x*z + y*z
    assert symmetric_poly(3, x, y, z) == x*y*z

def test_random_poly():
    poly = random_poly(x, 10, -100, 100, polys=False)

    assert Poly(poly).degree() == 10
    assert all(-100 <= coeff <= 100 for coeff in Poly(poly).coeffs()) is True

    poly = random_poly(x, 10, -100, 100, polys=True)

    assert poly.degree() == 10
    assert all(-100 <= coeff <= 100 for coeff in poly.coeffs()) is True

def test_fateman_poly_F_1():
    f,g,h = fateman_poly_F_1(1)
    F,G,H = dmp_fateman_poly_F_1(1, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

    f,g,h = fateman_poly_F_1(3)
    F,G,H = dmp_fateman_poly_F_1(3, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

def test_fateman_poly_F_2():
    f,g,h = fateman_poly_F_2(1)
    F,G,H = dmp_fateman_poly_F_2(1, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

    f,g,h = fateman_poly_F_2(3)
    F,G,H = dmp_fateman_poly_F_2(3, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

def test_fateman_poly_F_3():
    f,g,h = fateman_poly_F_3(1)
    F,G,H = dmp_fateman_poly_F_3(1, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

    f,g,h = fateman_poly_F_3(3)
    F,G,H = dmp_fateman_poly_F_3(3, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

