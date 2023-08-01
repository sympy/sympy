from hypothesis import given
from hypothesis import strategies as st
from sympy.abc import x
from sympy.polys.polytools import Poly


@given(
    coefficients1=st.lists(st.integers()),
    coefficients2=st.lists(st.integers()),
    coefficients3=st.lists(st.integers()),
)
def test_gcd_hypothesis(coefficients1, coefficients2, coefficients3):
    f = Poly(coefficients1, x, domain="ZZ")
    g = Poly(coefficients2, x, domain="ZZ")
    r = Poly(coefficients3, x, domain="ZZ")

    gcd_1 = f.gcd(g)
    gcd_2 = g.gcd(f)

    assert gcd_1 == gcd_2

    # multiply by r
    gcd_3 = g.gcd(f + r * g)

    assert gcd_1 == gcd_3


@given(
    coefficients1=st.lists(st.integers()),
    coefficients2=st.lists(st.integers()).filter(lambda x: any(x)),
)
def test_poly_hypothesis(coefficients1, coefficients2):
    # Integer case
    f_z = Poly(coefficients1, x, domain="ZZ")
    g_z = Poly(coefficients2, x, domain="ZZ")
    remainder_z = f_z.rem(g_z)
    assert g_z.degree() >= remainder_z.degree() or remainder_z.degree() == 0

    # Rational case
    f_q = Poly(coefficients1, x, domain="QQ")
    g_q = Poly(coefficients2, x, domain="QQ")
    remainder_q = f_q.rem(g_q)
    assert g_q.degree() >= remainder_q.degree() or remainder_q.degree() == 0
