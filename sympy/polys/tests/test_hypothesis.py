from hypothesis import given
from hypothesis import strategies as st
from sympy.testing.hypothesis import polynomial
from sympy import gcd, lcm



@given(
    f=polynomial(),
    g=polynomial(),
    r=polynomial(),
    x=st.integers(),
    y=st.integers(),
)
def test_gcd(f, g, r, x, y):
    gcd_1 = f.gcd(g)
    gcd_2 = g.gcd(f)

    assert gcd_1 == gcd_2

    # multiply by r
    gcd_3 = g.gcd(f + r * g)

    assert gcd_1 == gcd_3

    assert gcd(x, y) == gcd(y, x)
@given(
    f=polynomial(),
    g=polynomial(empty=False),
    h=polynomial(domain="QQ"),
    l=polynomial(empty=False, domain="QQ"),
)
def test_division(f, g, h, l):
    # Integer case
    z = f.rem(g)
    assert g.degree() >= z.degree() or z.degree() == 0

    # Rational case
    q = h.rem(l)
    assert l.degree() >= q.degree() or q.degree() == 0


@given(
    f=polynomial(),
    g=polynomial(),
)
def test_multiplication(f, g):
    h = f * g
    assert h.degree() == f.degree() + g.degree()
    assert h.LC() == f.LC() * g.LC()


@given(
    f=polynomial(),
    g=polynomial(),
)
def test_addition(f, g):
    h = f + g
    if h.degree() != -float("inf"):
        assert h.degree() == max(f.degree(), g.degree())


@given(
    f=polynomial(),
    g=polynomial(empty=False),
    x=st.integers(),
    y=st.integers(),
)
def test_lcm(f, g, x, y):
    assert f.lcm(g) == g.lcm(f)
    assert f * g == f.lcm(g) * f.gcd(g)
    if f.coeffs()[0] == 0 or g.coeffs()[0] == 0:
        assert f.lcm(g) == 0
    assert lcm(x, y) == lcm(y, x)


@given(
    f=polynomial(),
    g=polynomial(size=1),
)
def test_dispersion(f, g):
    assert f.dispersion() == f.dispersion(f)
    assert f.dispersion(g) == 0
