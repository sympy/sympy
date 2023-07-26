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
    # instantiate polynomials
    f = Poly(coefficients1, x, domain="ZZ")
    g = Poly(coefficients2, x, domain="ZZ")
    r = Poly(coefficients3, x, domain="ZZ")

    # get gdc
    gcd_1 = f.gcd(g)
    gcd_2 = g.gcd(f)

    assert gcd_1 == gcd_2

    # multiply by r
    gcd_3 = g.gcd(f + r * g)

    assert gcd_1 == gcd_3
