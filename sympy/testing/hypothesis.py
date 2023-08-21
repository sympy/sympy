from hypothesis import strategies as st
from hypothesis import assume
from sympy.abc import x
from sympy.polys.domains import ZZ, QQ
from sympy.polys.polytools import Poly


@st.composite
def coefficients(draw: st.DrawFn, empty=True, size=None):
    min_size = 0 if empty else 1
    # converting to the underlying type (ZZ).
    raw_l = draw(st.lists(st.integers(), min_size=min_size, max_size=size))
    l = [ZZ(i) for i in raw_l]
    # ensuring there is no leading zero or fully zero list.
    if len(l) > 0:
        assume(any(l))
        assume(l[0] != 0)
    return l


@st.composite
def coefficients_rational(draw: st.DrawFn, empty=True, size=None):
    min_size = 0 if empty else 1
    # converting to the underlying type (QQ).
    raw_l = draw(st.lists(st.fractions(), min_size=min_size, max_size=size))
    l = [QQ(i) for i in raw_l]
    # ensuring there is no leading zero or fully zero list.
    if len(l) > 0:
        assume(any(l))
        assume(l[0] != 0)
    return l


@st.composite
def polynomial(draw: st.DrawFn, empty=True, sparse=False, domain="ZZ", size=None):
    if not sparse:
        return Poly(draw(coefficients(empty, size=size)), x, domain=domain)
    # sparse polynomials are not created using coefficients.
    return []
