from hypothesis import given
from sympy.testing.hypothesis import polys


@given(
    f=polys(),
    g=polys(),
    r=polys(),
)
def test_gcd(f, g, r):
    gcd_1 = f.gcd(g)
    gcd_2 = g.gcd(f)

    assert gcd_1 == gcd_2

    # multiply by r
    gcd_3 = g.gcd(f + r * g)

    assert gcd_1 == gcd_3


@given(
    f=polys(),
    g=polys(empty=False),
    h=polys(domain="QQ"),
    l=polys(empty=False, domain="QQ"),
)
def test_division(f, g, h, l):
    # Integer case
    z = f.rem(g)
    assert g.degree() >= z.degree() or z.degree() == 0

    # Rational case
    q = h.rem(l)
    assert l.degree() >= q.degree() or q.degree() == 0


@given(
    f=polys(),
    g=polys(),
)
def test_multiplication(f, g):
    h = f * g
    assert h.degree() == f.degree() + g.degree()
    assert h.LC() == f.LC() * g.LC()


@given(
    f=polys(),
    g=polys(),
)
def test_addition(f, g):
    h = f + g
    if h.degree() != -float("inf"):
        assert h.degree() == max(f.degree(), g.degree())


@given(
    f=polys(),
    g=polys(empty=False),
)
def test_lcm(f, g):
    assert f.lcm(g) == g.lcm(f)
    assert f * g == f.lcm(g) * f.gcd(g)
    if f.coeffs()[0] == 0 or g.coeffs()[0] == 0:
        assert f.lcm(g) == 0


@given(
    f=polys(),
    g=polys(degree=1),
)
def test_dispersion(f, g):
    assert f.dispersion() == f.dispersion(f)
    assert f.dispersion(g) == 0
