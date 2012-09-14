"""Tests for tools and arithmetics for monomials of distributed polynomials. """

from sympy.polys.monomialtools import (
    monomials, monomial_count,
    monomial_key, lex, grlex, grevlex, ilex, igrlex, igrevlex,
    monomial_mul, monomial_div,
    monomial_gcd, monomial_lcm,
    monomial_max, monomial_min,
    monomial_divides,
    Monomial,
    InverseOrder, ProductOrder, build_product_order,
    LexOrder
)

from sympy.polys.polyerrors import ExactQuotientFailed

from sympy.abc import a, b, c, x, y, z
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

    assert lex.is_global is True
    assert lex == LexOrder()
    assert lex != grlex

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

    assert grlex.is_global is True

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

    assert grevlex.is_global is True

def test_InverseOrder():
    ilex = InverseOrder(lex)
    igrlex = InverseOrder(grlex)

    assert ilex((1,2,3)) > ilex((2, 0, 3))
    assert igrlex((1, 2, 3)) < igrlex((0, 2, 3))
    assert str(ilex) == "ilex"
    assert str(igrlex) == "igrlex"
    assert ilex.is_global is False
    assert igrlex.is_global is False
    assert ilex != igrlex
    assert ilex == InverseOrder(LexOrder())

def test_ProductOrder():
    P = ProductOrder((grlex, lambda m: m[:2]), (grlex, lambda m: m[2:]))
    assert P((1, 3, 3, 4, 5)) > P((2, 1, 5, 5, 5))
    assert str(P) == "ProductOrder(grlex, grlex)"
    assert P.is_global is True
    assert ProductOrder((grlex, None), (ilex, None)).is_global is None
    assert ProductOrder((igrlex, None), (ilex, None)).is_global is False

def test_monomial_key():
    assert monomial_key() == lex

    assert monomial_key('lex') == lex
    assert monomial_key('grlex') == grlex
    assert monomial_key('grevlex') == grevlex

    raises(ValueError, lambda: monomial_key('foo'))
    raises(ValueError, lambda: monomial_key(1))

def test_build_product_order():
    from sympy.abc import x, y, z, t
    assert build_product_order((("grlex", x, y), ("grlex", z, t)), [x, y, z, t]) \
        ((4, 5, 6, 7)) == ((9, (4, 5)), (13, (6, 7)))

    assert build_product_order((("grlex", x, y), ("grlex", z, t)), [x, y, z, t]) == \
               build_product_order((("grlex", x, y), ("grlex", z, t)), [x, y, z, t])
    assert (build_product_order((("grlex", x, y), ("grlex", z, t)), [x, y, z, t]) != \
               build_product_order((("grlex", x, y), ("grlex", z, t)), [x, y, z, t])) \
           is False

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

def test_monomial_divides():
    assert monomial_divides((1,2,3), (4,5,6)) is True
    assert monomial_divides((1,2,3), (0,5,6)) is False

def test_Monomial():
    m = Monomial((3, 4, 1), (x, y, z))
    n = Monomial((1, 2, 0), (x, y, z))

    assert m.as_expr() == x**3*y**4*z
    assert n.as_expr() == x**1*y**2

    assert m.as_expr(a, b, c) == a**3*b**4*c
    assert n.as_expr(a, b, c) == a**1*b**2

    assert m.exponents == (3, 4, 1)
    assert m.gens == (x, y, z)

    assert n.exponents == (1, 2, 0)
    assert n.gens == (x, y, z)

    assert m == (3, 4, 1)
    assert n != (3, 4, 1)
    assert m != (1, 2, 0)
    assert n == (1, 2, 0)

    assert m[0] == m[-3] == 3
    assert m[1] == m[-2] == 4
    assert m[2] == m[-1] == 1

    assert n[0] == n[-3] == 1
    assert n[1] == n[-2] == 2
    assert n[2] == n[-1] == 0

    assert m[:2] == (3, 4)
    assert n[:2] == (1, 2)

    assert m*n == Monomial((4, 6, 1))
    assert m/n == Monomial((2, 2, 1))

    assert m*(1, 2, 0) == Monomial((4, 6, 1))
    assert m/(1, 2, 0) == Monomial((2, 2, 1))

    assert m.gcd(n) == Monomial((1, 2, 0))
    assert m.lcm(n) == Monomial((3, 4, 1))

    assert m.gcd((1, 2, 0)) == Monomial((1, 2, 0))
    assert m.lcm((1, 2, 0)) == Monomial((3, 4, 1))

    assert m**0 == Monomial((0, 0, 0))
    assert m**1 == m
    assert m**2 == Monomial((6, 8, 2))
    assert m**3 == Monomial((9,12, 3))

    raises(ExactQuotientFailed, lambda: m/Monomial((5, 2, 0)))
