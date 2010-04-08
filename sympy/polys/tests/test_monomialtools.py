"""Tests for tools and arithmetics for monomials of distributed polynomials. """

from sympy.polys.monomialtools import (
    monomials, monomial_count,
    monomial_lex_key, monomial_grlex_key, monomial_grevlex_key, monomial_key,
    monomial_lex_cmp, monomial_grlex_cmp, monomial_grevlex_cmp, monomial_cmp,
    monomial_mul, monomial_div,
    monomial_gcd, monomial_lcm,
    monomial_max, monomial_min,
)

from sympy.abc import x, y
from sympy.utilities.pytest import raises

def test_monomials():
    assert sorted(monomials([], 0)) == [1]
    assert sorted(monomials([], 1)) == [1]
    assert sorted(monomials([], 2)) == [1]
    assert sorted(monomials([], 3)) == [1]

    assert sorted(monomials([x], 0)) == [1]
    assert sorted(monomials([x], 1)) == [1, x]
    assert sorted(monomials([x], 2)) == [1, x, x**2]
    assert sorted(monomials([x], 3)) == [1, x, x**2, x**3]

    assert sorted(monomials([x, y], 0)) == [1]
    assert sorted(monomials([x, y], 1)) == [1, x, y]
    assert sorted(monomials([x, y], 2)) == [1, x, y, x**2, y**2, x*y]
    assert sorted(monomials([x, y], 3)) == [1, x, y, x**2, x**3, y**2, y**3, x*y, x*y**2, y*x**2]

def test_monomial_count():
    assert monomial_count(2, 2) == 6
    assert monomial_count(2, 3) == 10

def test_monomial_lex_key():
    assert monomial_lex_key((1,2,3)) == (1,2,3)

def test_monomial_grlex_key():
    assert monomial_grlex_key((1,2,3)) == (6, (1,2,3))

def test_monomial_grevlex_key():
    assert monomial_grevlex_key((1,2,3)) == (6, (3,2,1))

def test_monomial_key():
    assert monomial_key('lex') == monomial_lex_key
    assert monomial_key('grlex') == monomial_grlex_key
    assert monomial_key('grevlex') == monomial_grevlex_key

    raises(ValueError, "monomial_key('foo')")
    raises(ValueError, "monomial_key(1)")

def test_monomial_lex_cmp():
    assert monomial_lex_cmp((1,2,3), (1,2,3)) == 0

    assert monomial_lex_cmp((2,2,3), (1,2,3)) == 1
    assert monomial_lex_cmp((1,3,3), (1,2,3)) == 1
    assert monomial_lex_cmp((1,2,4), (1,2,3)) == 1

    assert monomial_lex_cmp((0,2,3), (1,2,3)) == -1
    assert monomial_lex_cmp((1,1,3), (1,2,3)) == -1
    assert monomial_lex_cmp((1,2,2), (1,2,3)) == -1

def test_monomial_grlex_cmp():
    assert monomial_grlex_cmp((1,2,3), (1,2,3)) == 0

    assert monomial_grlex_cmp((2,2,3), (1,2,3)) == 1
    assert monomial_grlex_cmp((1,3,3), (1,2,3)) == 1
    assert monomial_grlex_cmp((1,2,4), (1,2,3)) == 1

    assert monomial_grlex_cmp((0,2,3), (1,2,3)) == -1
    assert monomial_grlex_cmp((1,1,3), (1,2,3)) == -1
    assert monomial_grlex_cmp((1,2,2), (1,2,3)) == -1

    assert monomial_grlex_cmp((2,2,3), (1,2,4)) == 1
    assert monomial_grlex_cmp((1,3,3), (1,2,4)) == 1

    assert monomial_grlex_cmp((0,2,3), (1,2,2)) == -1
    assert monomial_grlex_cmp((1,1,3), (1,2,2)) == -1

def test_monomial_grevlex_cmp():
    assert monomial_grevlex_cmp((1,2,3), (1,2,3)) == 0

    assert monomial_grevlex_cmp((2,2,3), (1,2,3)) == 1
    assert monomial_grevlex_cmp((1,3,3), (1,2,3)) == 1
    assert monomial_grevlex_cmp((1,2,4), (1,2,3)) == 1

    assert monomial_grevlex_cmp((0,2,3), (1,2,3)) == -1
    assert monomial_grevlex_cmp((1,1,3), (1,2,3)) == -1
    assert monomial_grevlex_cmp((1,2,2), (1,2,3)) == -1

    assert monomial_grevlex_cmp((2,2,3), (1,2,4)) == 1
    assert monomial_grevlex_cmp((1,3,3), (1,2,4)) == 1

    assert monomial_grevlex_cmp((0,2,3), (1,2,2)) == -1
    assert monomial_grevlex_cmp((1,1,3), (1,2,2)) == -1

def test_monomial_cmp():
    assert monomial_cmp('lex') == monomial_lex_cmp
    assert monomial_cmp('grlex') == monomial_grlex_cmp
    assert monomial_cmp('grevlex') == monomial_grevlex_cmp

    raises(ValueError, "monomial_cmp('unknown')")

def test_monomial_mul():
    assert monomial_mul((3,4,1), (1,2,0)) == (4,6,1)

def test_monomial_div():
    assert monomial_div((3,4,1), (1,2,0)) == (2,2,1)

def test_monomial_gcd():
    assert monomial_gcd((3,4,1), (1,2,0)) == (1,2,0)

def test_monomial_lcm():
    assert monomial_lcm((3,4,1), (1,2,0)) == (3,4,1)

def test_monomial_max():
    assert monomial_max((3,4,5), (0,5,1), (6,3,9)) == (6,5,9)

def test_monomial_min():
    assert monomial_min((3,4,5), (0,5,1), (6,3,9)) == (0,3,1)

