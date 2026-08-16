from __future__ import annotations
from itertools import product
from hypothesis import assume, given, reject
from hypothesis import settings, strategies as st
from hypothesis.strategies import composite

from sympy.abc import x
from sympy import ZZ, QQ
from sympy.polys.heuristicgcd import heugcd
from sympy.polys.polyerrors import HeuristicGCDFailed
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
from sympy.polys.rings import ring
from sympy.polys.sparsetools import smp_mul
from sympy.polys.zippel import smp_zippel_gcd
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


@settings(max_examples=25)
@given(
    f=dup_any(dom=QQ, max_len=20),
    g=dup_zero_TC(dom=QQ, max_len=20),
    n=st.integers(min_value=3, max_value=20),
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


def test_dup_series_mul_rational_karatsuba_boundary():
    # dup_series_mul switches to Karatsuba once min(len(f), len(g)) >= 100.
    # Keep this as a separate test so the compose Hypothesis test can stay small
    # and avoid paying for large exact dup_compose() oracles.
    f = dup_from_list([QQ(1)] + [QQ(1, 2)] * 98 + [QQ(1, 3)], QQ)
    g = dup_from_list([QQ(2)] + [QQ(-1, 3)] * 98 + [QQ(1, 5)], QQ)

    base = _dup_series_mul_base(f, g, 100, QQ)
    karatsuba = _dup_series_mul_karatsuba(f, g, 100, QQ)
    series_mul = dup_series_mul(f, g, 100, QQ)

    assert series_mul == karatsuba
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


nonzero_coefficients = st.one_of(
    st.integers(min_value=-1000, max_value=-1),
    st.integers(min_value=1, max_value=1000),
)


@composite
def smp_polynomials(draw, n, max_degree, nonconstant=False):
    monomials = [
        mon
        for mon in product(range(max_degree + 1), repeat=n)
        if sum(mon) <= max_degree
    ]

    max_terms = min(10, len(monomials))
    polynomial = draw(
        st.dictionaries(
            keys=st.sampled_from(monomials),
            values=nonzero_coefficients,
            min_size=1,
            max_size=max_terms,
        )
    )

    if nonconstant:
        nonconstant_monomials = [mon for mon in monomials if any(mon)]
        monomial = draw(st.sampled_from(nonconstant_monomials))
        coefficient = draw(nonzero_coefficients)
        polynomial[monomial] = coefficient

    return {monomial: ZZ(coefficient) for monomial, coefficient in polynomial.items()}


@composite
def zippel_cases(draw, min_n):

    n = draw(st.integers(min_value=min_n, max_value=5))
    max_degree = draw(st.integers(min_value=1, max_value=10))

    a = draw(smp_polynomials(n, max_degree))
    b = draw(smp_polynomials(n, max_degree))
    common = draw(
        smp_polynomials(
            n,
            max_degree,
            nonconstant=True,
        )
    )

    return n, max_degree, a, b, common


def _coefficients_wrt_main(poly):
    coefficients = {}

    for monomial, coefficient in poly.items():
        degree = monomial[0]
        coefficient_poly = coefficients.setdefault(degree, {})
        coefficient_poly[(0,) + monomial[1:]] = coefficient

    return list(coefficients.values())


def _content_gcd_wrt_main(f, g, ring):
    coefficients = _coefficients_wrt_main(f) + _coefficients_wrt_main(g)

    content = ring.from_dict(coefficients[0])

    for coefficient in coefficients[1:]:
        content = heugcd(
            content,
            ring.from_dict(coefficient),
        )[0]

        if content.is_ground:
            break

    return content


def _normalize_gcd_result(result):
    gcd, cff, cfg = result

    if gcd.LC < 0:
        gcd = -gcd
        cff = -cff
        cfg = -cfg

    return gcd, cff, cfg


@given(
    case=zippel_cases(1),
    random_state=st.random_module(),
)
def test_smp_zippel_gcd(case, random_state):
    n, _, a, b, common = case

    R, *_ = ring(f"x:{n}", ZZ)

    f = smp_mul(a, common, ZZ, n)
    g = smp_mul(b, common, ZZ, n)

    try:
        content_gcd = _content_gcd_wrt_main(f, g, R)
    except HeuristicGCDFailed:
        reject()

    # Zippel currently expects no common polynomial content
    # with respect to the main variable.
    assume(content_gcd.is_ground)

    try:
        expected = heugcd(R(f), R(g))
    except HeuristicGCDFailed:
        reject()

    gcd, cff, cfg = smp_zippel_gcd(f, g, n)

    expected = _normalize_gcd_result(expected)

    result = _normalize_gcd_result((R(gcd), R(cff), R(cfg)))

    assert result == expected


@given(
    case=zippel_cases(1),
    random_state=st.random_module(),
)
def test_smp_zippel_gcd_independent(case, random_state):
    n, _, f, g, _ = case
    R, *_ = ring(f"x:{n}", ZZ)

    try:
        content_gcd = _content_gcd_wrt_main(f, g, R)
    except HeuristicGCDFailed:
        reject()

    assume(content_gcd.is_ground)

    try:
        expected = heugcd(R.from_dict(f), R.from_dict(g))
    except HeuristicGCDFailed:
        reject()

    gcd, cff, cfg = smp_zippel_gcd(f, g, n)

    expected = _normalize_gcd_result(expected)
    result = _normalize_gcd_result((R(gcd), R(cff), R(cfg)))

    assert result == expected


@given(
    case=zippel_cases(2),
    random_state=st.random_module(),
)
def test_smp_zippel_gcd_monic(case, random_state):
    n, max_degree, a, _, h = case
    R, *_ = ring(f"x:{n}", ZZ)

    # building cofactor b, coprime to a
    zm = (0,) * n
    b = a.copy()
    b[zm] = b.get(zm, ZZ.zero) + ZZ.one
    if not b[zm]:
        del b[zm]

    # building monic gcd
    d = max_degree + 1
    h[(d,) + (0,) * (n - 1)] = ZZ.one

    f = smp_mul(a, h, ZZ, n)
    g = smp_mul(b, h, ZZ, n)
    gcd, _, _ = smp_zippel_gcd(f, g, n)
    assert gcd == h


@given(
    case=zippel_cases(3),
    random_state=st.random_module(),
)
def test_smp_zippel_gcd_pseudomonic(case, random_state):
    n, max_degree, a, _, h = case
    R, *_ = ring(f"x:{n}", ZZ)

    # building cofactor b, coprime to a
    zm = (0,) * n
    b = a.copy()
    b[zm] = b.get(zm, ZZ.zero) + ZZ.one
    if not b[zm]:
        del b[zm]

    # building pseudomonic gcd
    d = max_degree + 1
    h[(d, 1) + (0,) * (n - 2)] = ZZ.one
    h[(max_degree,) + (0,) * (n - 1)] = ZZ.one

    f = smp_mul(a, h, ZZ, n)
    g = smp_mul(b, h, ZZ, n)
    gcd, _, _ = smp_zippel_gcd(f, g, n)
    assert gcd == h


@given(
    case=zippel_cases(3),
    random_state=st.random_module(),
)
def test_smp_zippel_gcd_not_monic(case, random_state):
    n, max_degree, a, _, h = case
    R, *_ = ring(f"x:{n}", ZZ)

    # building cofactor b, coprime to a
    zm = (0,) * n
    b = a.copy()
    b[zm] = b.get(zm, ZZ.zero) + ZZ.one
    if not b[zm]:
        del b[zm]

    # building not monic gcd
    d = max_degree + 1
    h[(d,) + (0,) * (n - 1)] = ZZ.one
    h[(d, 1) + (0,) * (n - 2)] = ZZ.one
    h[(max_degree,) + (0,) * (n - 1)] = ZZ.one

    f = smp_mul(a, h, ZZ, n)
    g = smp_mul(b, h, ZZ, n)
    gcd, _, _ = smp_zippel_gcd(f, g, n)
    assert gcd == h
