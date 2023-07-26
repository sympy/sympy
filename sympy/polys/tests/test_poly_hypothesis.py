from hypothesis import given
from hypothesis import assume
from hypothesis import strategies as st
from sympy.abc import x
from sympy.polys.polytools import Poly


@given(coefficients1=st.lists(st.integers()), coefficients2=st.lists(st.integers()))
def test_poly_hypothesis(coefficients1, coefficients2):
    # non-zero division check.
    assume(len(coefficients2) * [0] != coefficients2)

    # Integer case
    f_z = Poly(coefficients1, x, domain="ZZ")
    g_z = Poly(coefficients2, x, domain="ZZ")
    remainder_z = f_z.rem(g_z)
    assert g_z.degree() > remainder_z.degree() or remainder_z.degree() == 0

    # Rational case
    f_q = Poly(coefficients1, x, domain="QQ")
    g_q = Poly(coefficients2, x, domain="QQ")
    remainder_q = f_q.rem(g_q)
    assert g_q.degree() > remainder_q.degree() or remainder_q.degree() == 0
