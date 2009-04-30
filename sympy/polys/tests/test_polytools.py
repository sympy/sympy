"""Tests for user-friendly public interface to polynomial functions. """

from sympy.polys.polytools import (
    Poly,
    pdiv, prem, pquo, pexquo,
    div, rem, quo, exquo,
    half_gcdex, gcdex, invert,
    subresultants,
    resultant, discriminant,
    cofactors, gcd, lcm,
    monic, content, primitive,
    compose, decompose,
    sqf_part, sqf_list, sqf,
    factor_list, factor,
    cancel, sturm,
    groebner,
)

from sympy import S, symbols, sqrt, exp

x, y = symbols('x,y')

def test__init_poly_from_poly():
    assert Poly(x*y, x, y).gens == (x, y)
    assert Poly(x*y, y, x).gens == (y, x)

    assert Poly(Poly(x*y, x, y), y, x).gens == (y, x)

def test__init_poly_from_basic():
    assert Poly(0).is_Poly == False
    assert Poly(1).is_Poly == False

    assert Poly(0, x).is_Poly == True
    assert Poly(1, x).is_Poly == True

def test_pdiv():
    pass

def test_prem():
    pass

def test_pquo():
    pass

def test_pexquo():
    pass

def test_div():
    pass

def test_rem():
    pass

def test_quo():
    pass

def test_exquo():
    pass

def test_half_gcdex():
    pass

def test_gcdex():
    pass

def test_invert():
    pass

def test_subresultants():
    pass

def test_resultant():
    pass

def test_discriminant():
    pass

def test_cofactors():
    pass

def test_gcd():
    pass

def test_lcm():
    pass

def test_monic():
    pass

def test_content():
    pass

def test_primitive():
    pass

def test_compose():
    pass

def test_decompose():
    pass

def test_sqf_part():
    pass

def test_sqf_list():
    pass

def test_sqf():
    pass

def test_factor_list():
    pass

def test_factor():
    pass

def test_sturm():
    pass

def test_cancel():
    assert cancel(0) == 0
    assert cancel(7) == 7
    assert cancel(x) == x

    f = (x**2 - 2)/(x + sqrt(2))

    assert cancel(f) == f
    assert cancel(f, greedy=False) == x - sqrt(2)

    f = (x**2 - 2)/(x - sqrt(2))

    assert cancel(f) == f
    assert cancel(f, greedy=False) == x + sqrt(2)

    assert cancel((x**2-y)/(x-y)) == 1/(x - y)*(x**2 - y)

    assert cancel((x**2-y**2)/(x-y), x) == x + y
    assert cancel((x**2-y**2)/(x-y), y) == x + y
    assert cancel((x**2-y**2)/(x-y)) == x + y

    assert cancel((x**3-1)/(x**2-1)) == (x**2+x+1)/(x+1)
    assert cancel((x**3/2-S(1)/2)/(x**2-1)) == (x**2+x+1)/(2*x+2)

    assert cancel((exp(2*x) + 2*exp(x) + 1)/(exp(x) + 1)) == exp(x) + 1

def test_groebner():
    pass

