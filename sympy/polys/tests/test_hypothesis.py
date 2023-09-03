from hypothesis import given
from sympy.testing.hypothesis import polys, lattice_axioms_singular, lattice_axioms_dual
from sympy.polys.polytools import lcm, gcd


@given(
    f=polys(),
    g=polys(),
    h=polys(),
)
def test_gcd(f, g, h):
    assert lattice_axioms_singular([f, g, h], gcd)

    # multiply by h
    gcd_h = g.gcd(f + h * g)
    assert f.gcd(g) == gcd_h


@given(f=polys(), g=polys(), h=polys())
def test_lcm(f, g, h):
    assert lattice_axioms_singular([f, g, h], lcm)
    assert f * g == f.lcm(g) * f.gcd(g)


@given(f=polys(), g=polys(), h=polys())
def test_lcm_gcd(f, g, h):
    assert lattice_axioms_dual([f, g, h], [lcm, gcd])


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
    g=polys(degree=1),
)
def test_dispersion(f, g):
    assert f.dispersion() == f.dispersion(f)
    assert f.dispersion(g) == 0
