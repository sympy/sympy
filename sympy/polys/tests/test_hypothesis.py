from hypothesis import given
from hypothesis import settings, strategies as st
from hypothesis.strategies import composite

from sympy.abc import x
from sympy import ZZ, QQ
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
def dup_unit(draw, dom=ZZ, min_len=3, max_len=200):
    if dom == ZZ:
        elem_strategy = st.integers()
    elif dom == QQ:
        elem_strategy = st.tuples(
            st.integers(min_value=-10, max_value=10),
            st.integers(min_value=1, max_value=10),
        ).map(lambda t: QQ(t[0], t[1]))

    lst = draw(st.lists(elem_strategy, min_size=min_len, max_size=max_len))

    lst[-1] = dom(0)
    if dom == ZZ:
        lst[-2] = dom(1)
    elif dom == QQ:
        if dom.is_zero(lst[-2]):
            lst[-2] = dom(1)

    return lst


@composite
def dup_any(draw, dom=ZZ, min_len=3, max_len=200):
    if dom == ZZ:
        elem_strategy = st.integers()
    elif dom == QQ:
        elem_strategy = st.tuples(st.integers(), st.integers(min_value=1)).map(
            lambda t: QQ(t[0], t[1])
        )

    return draw(st.lists(elem_strategy, min_size=min_len, max_size=max_len))


@composite
def dup_zero_TC(draw, dom=ZZ, min_len=3, max_len=200):
    if dom == ZZ:
        elem_strategy = st.integers()
    elif dom == QQ:
        elem_strategy = st.tuples(st.integers(), st.integers(min_value=1)).map(
            lambda t: QQ(t[0], t[1])
        )

    lst = draw(st.lists(elem_strategy, min_size=min_len, max_size=max_len))

    lst[-1] = dom(0)

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


@given(f=dup_any(), g=dup_zero_TC(), n=st.integers(min_value=3, max_value=200))
def test_dup_series_compose_int(f, g, n):
    expected = dup_truncate(dup_compose(f, g, ZZ), n, ZZ)
    assert dup_series_compose(f, g, n, ZZ) == expected


@slow
@given(f=dup_unit(), n=st.integers(min_value=3, max_value=200))
def test_dup_series_reversion_int(f, n):
    rev = dup_series_reversion(f, n, ZZ)
    comp = dup_series_compose(rev, f, n, ZZ)
    expected = dup_from_list([1, 0], ZZ)
    assert comp == expected


@given(f=dup_any(), g=dup_any(), n=st.integers(min_value=3, max_value=200))
def test_dup_series_mul_int(f, g, n):
    base = _dup_series_mul_base(f, g, n, ZZ)
    karatsuba = _dup_series_mul_karatsuba(f, g, n, ZZ)
    series_mul = dup_series_mul(f, g, n, ZZ)
    Dup_mul = dup_truncate(dup_mul(f, g, ZZ), n, ZZ)

    assert series_mul == Dup_mul
    assert base == karatsuba


@given(
    f=dup_any(),
    pow=st.integers(min_value=0, max_value=10),
    n=st.integers(min_value=3, max_value=200),
)
def test_dup_series_pow_int(f, pow, n):
    Dup_pow = dup_truncate(dup_pow(f, pow, ZZ), n, ZZ)
    series_pow = dup_series_pow(f, pow, n, ZZ)

    assert series_pow == Dup_pow


@given(
    f=dup_any(dom=QQ), g=dup_zero_TC(dom=QQ), n=st.integers(min_value=3, max_value=200)
)
def test_dup_series_compose_rational(f, g, n):
    expected = dup_truncate(dup_compose(f, g, QQ), n, QQ)
    assert dup_series_compose(f, g, n, QQ) == expected


@slow
@given(f=dup_unit(dom=QQ), n=st.integers(min_value=3, max_value=150))
@settings(max_examples=100)
def test_dup_series_reversion_rational(f, n):
    rev = dup_series_reversion(f, n, QQ)
    comp = dup_series_compose(rev, f, n, QQ)
    expected = dup_from_list([1, 0], QQ)
    assert comp == expected


@given(f=dup_any(dom=QQ), g=dup_any(dom=QQ), n=st.integers(min_value=3, max_value=100))
def test_dup_series_mul_rational(f, g, n):
    base = _dup_series_mul_base(f, g, n, QQ)
    karatsuba = _dup_series_mul_karatsuba(f, g, n, QQ)
    series_mul = dup_series_mul(f, g, n, QQ)
    Dup_mul = dup_truncate(dup_mul(f, g, QQ), n, QQ)

    assert series_mul == Dup_mul
    assert base == karatsuba


@given(
    f=dup_any(dom=QQ),
    pow=st.integers(min_value=0, max_value=10),
    n=st.integers(min_value=3, max_value=200),
)
def test_dup_series_pow_rational(f, pow, n):
    Dup_pow = dup_truncate(dup_pow(f, pow, QQ), n, QQ)
    series_pow = dup_series_pow(f, pow, n, QQ)

    assert series_pow == Dup_pow
