from hypothesis import given
from hypothesis import strategies as st
from hypothesis.strategies import composite

from sympy.abc import x
from sympy import ZZ
from sympy.polys.polytools import Poly
from sympy.polys.densebasic import dup_truncate, dup_from_list
from sympy.polys.densetools import dup_compose, dup_series_compose, dup_series_reversion
from sympy.polys.densearith import (
    dup_mul,
    dup_series_mul,
    _dup_series_mul_base,
    _dup_series_mul_karatsuba,
    dup_pow,
    dup_series_pow,
)
from sympy.testing.pytest import slow


def polys(*, nonzero=False, domain="ZZ"):
    # This is a simple strategy, but sufficient the tests below
    elems = {"ZZ": st.integers(), "QQ": st.fractions()}
    coeff_st = st.lists(elems[domain])
    if nonzero:
        coeff_st = coeff_st.filter(any)
    return st.builds(Poly, coeff_st, st.just(x), domain=st.just(domain))


@composite
def dup(draw, unit=False, min_len=3, max_len=200):
    lst = draw(st.lists(st.integers(), min_size=min_len, max_size=max_len))

    if unit:
        lst[-1] = ZZ(0)
        lst[-2] = ZZ(1)

    return lst


@given(f=polys(), g=polys(), r=polys())
def test_gcd_hypothesis(f, g, r):
    gcd_1 = f.gcd(g)
    gcd_2 = g.gcd(f)
    assert gcd_1 == gcd_2

    # multiply by r
    gcd_3 = g.gcd(f + r * g)
    assert gcd_1 == gcd_3


@given(f_z=polys(), g_z=polys(nonzero=True))
def test_poly_hypothesis_integers(f_z, g_z):
    remainder_z = f_z.rem(g_z)
    assert g_z.degree() >= remainder_z.degree() or remainder_z.degree() == 0


@given(f_q=polys(domain="QQ"), g_q=polys(nonzero=True, domain="QQ"))
def test_poly_hypothesis_rationals(f_q, g_q):
    remainder_q = f_q.rem(g_q)
    assert g_q.degree() >= remainder_q.degree() or remainder_q.degree() == 0


@given(f=dup(unit=True), g=dup(unit=True), n=st.integers(min_value=3, max_value=200))
def test_dup_series_compose(f, g, n):
    expected = dup_truncate(dup_compose(f, g, ZZ), n, ZZ)
    assert dup_series_compose(f, g, n, ZZ) == expected


@slow
@given(f=dup(unit=True), n=st.integers(min_value=3, max_value=200))
def test_dup_series_reversion(f, n):
    rev = dup_series_reversion(f, n, ZZ)
    comp = dup_series_compose(rev, f, n, ZZ)
    expected = dup_from_list([1, 0], ZZ)
    assert comp == expected


@given(f=dup(), g=dup(), n=st.integers(min_value=3, max_value=200))
def test_dup_series_mul(f, g, n):
    base = _dup_series_mul_base(f, g, n, ZZ)
    karatsuba = _dup_series_mul_karatsuba(f, g, n, ZZ)
    series_mul = dup_series_mul(f, g, n, ZZ)
    Dup_mul = dup_truncate(dup_mul(f, g, ZZ), n, ZZ)

    assert series_mul == Dup_mul
    assert base == karatsuba


@given(
    f=dup(),
    pow=st.integers(min_value=0, max_value=10),
    n=st.integers(min_value=3, max_value=200),
)
def test_dup_series_pow(f, pow, n):
    Dup_pow = dup_truncate(dup_pow(f, pow, ZZ), n, ZZ)
    series_pow = dup_series_pow(f, pow, n, ZZ)

    assert series_pow == Dup_pow
