"""Tests for tools and arithmetics for monomials of distributed polynomials. """

from sympy.polys.monomialtools import (
    monomials, monomial_count,
    monomial_key, lex, grlex, grevlex,
    monomial_mul, monomial_div,
    monomial_gcd, monomial_lcm,
    monomial_max, monomial_min,
    Monomial,
)

from sympy.polys.polyerrors import ExactQuotientFailed

from sympy.abc import x, y, z
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

def test_lex_order():
    assert lex((1,2,3)) == (1,2,3)
    assert str(lex) == 'lex'

    assert lex((1,2,3)) == lex((1,2,3))

    assert lex((2,2,3)) > lex((1,2,3))
    assert lex((1,3,3)) > lex((1,2,3))
    assert lex((1,2,4)) > lex((1,2,3))

    assert lex((0,2,3)) < lex((1,2,3))
    assert lex((1,1,3)) < lex((1,2,3))
    assert lex((1,2,2)) < lex((1,2,3))

def test_grlex_order():
    assert grlex((1,2,3)) == (6, (1,2,3))
    assert str(grlex) == 'grlex'

    assert grlex((1,2,3)) == grlex((1,2,3))

    assert grlex((2,2,3)) > grlex((1,2,3))
    assert grlex((1,3,3)) > grlex((1,2,3))
    assert grlex((1,2,4)) > grlex((1,2,3))

    assert grlex((0,2,3)) < grlex((1,2,3))
    assert grlex((1,1,3)) < grlex((1,2,3))
    assert grlex((1,2,2)) < grlex((1,2,3))

    assert grlex((2,2,3)) > grlex((1,2,4))
    assert grlex((1,3,3)) > grlex((1,2,4))

    assert grlex((0,2,3)) < grlex((1,2,2))
    assert grlex((1,1,3)) < grlex((1,2,2))

def test_grevlex_order():
    assert grevlex((1,2,3)) == (6, (-3,-2,-1))
    assert str(grevlex) == 'grevlex'

    assert grevlex((1,2,3)) == grevlex((1,2,3))

    assert grevlex((2,2,3)) > grevlex((1,2,3))
    assert grevlex((1,3,3)) > grevlex((1,2,3))
    assert grevlex((1,2,4)) > grevlex((1,2,3))

    assert grevlex((0,2,3)) < grevlex((1,2,3))
    assert grevlex((1,1,3)) < grevlex((1,2,3))
    assert grevlex((1,2,2)) < grevlex((1,2,3))

    assert grevlex((2,2,3)) > grevlex((1,2,4))
    assert grevlex((1,3,3)) > grevlex((1,2,4))

    assert grevlex((0,2,3)) < grevlex((1,2,2))
    assert grevlex((1,1,3)) < grevlex((1,2,2))

    assert grevlex((0,1,1)) > grevlex((0,0,2))
    assert grevlex((0,3,1)) < grevlex((2,2,1))

def test_monomial_key():
    assert monomial_key() == lex

    assert monomial_key('lex') == lex
    assert monomial_key('grlex') == grlex
    assert monomial_key('grevlex') == grevlex

    raises(ValueError, "monomial_key('foo')")
    raises(ValueError, "monomial_key(1)")

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

def test_Monomial():
    assert Monomial(1, 2, 3).data == (1, 2, 3)

    m, n = Monomial(3, 4, 1), Monomial(1, 2, 0)

    assert m*n == Monomial(4, 6, 1)
    assert m/n == Monomial(2, 2, 1)

    assert m.gcd(n) == Monomial(1, 2, 0)
    assert m.lcm(n) == Monomial(3, 4, 1)

    a, b, c = [ Monomial(*monom) for monom in [(3,4,5), (0,5,1), (6,3,9)] ]

    assert Monomial.max(a, b, c) == Monomial(6, 5, 9)
    assert Monomial.min(a, b, c) == Monomial(0, 3, 1)

    n = Monomial(5, 2, 0)

    raises(ExactQuotientFailed, "m / n")

    assert n.as_expr(x, y, z) == x**5*y**2
