from sympy import symbols, ring, ZZ, QQ, lex, ExactQuotientFailed
from sympy.polys.fields import FracField
from sympy.polys.laurent import laurent_ring, LaurentPolyRing, LaurentPolyElement

from sympy.testing.pytest import raises


def test_LaurentRing_init():
    x, y = symbols("x y")
    Rp, _, _ = ring("x y", ZZ)
    R1, xr, yr = laurent_ring("x y", ZZ)
    R2 = LaurentPolyRing([x, y], ZZ, "lex")
    R3, _, _ = laurent_ring([x, y], ZZ)

    assert R1 == R2 == R3

    for R in [R1, R2, R3]:
        assert R.numer_ring == Rp
        assert R.domain == ZZ
        assert R.gens == (xr, yr)
        assert R.symbols == (x, y)
        assert R.order == lex
        assert R.zero == LaurentPolyElement(R, Rp.zero, Rp.one)
        assert R.one == LaurentPolyElement(R, Rp.one, Rp.one)
        assert R.x == xr
        assert R.y == yr


def test_LaurentRing_eq_hash():
    x, y = symbols("x y")
    R1, _, _ = laurent_ring("x y", ZZ)
    R2, _, _ = laurent_ring("x y", ZZ)
    R3, _, _ = laurent_ring("y x", ZZ)
    R4, _, _ = laurent_ring("x y", QQ)
    R5, _, _ = laurent_ring("x y", ZZ, order='grevlex')
    Rs = [R1, R2, R3, R4, R5, ZZ, QQ, ZZ[x,y]]

    for i, Ri in enumerate(Rs):
        for j, Rj in enumerate(Rs):
            if i == j or {i, j} == {0, 1}:
                assert hash(Ri) == hash(Rj)
                assert (Ri == Rj) is True
                assert (Ri != Rj) is False
            else:
                assert (Ri == Rj) is False
                assert (Ri != Rj) is True

    assert (R1 == QQ) is False
    assert (R1 != QQ) is True


def test_LaurentRing_str():
    R, x, y = laurent_ring("x y", ZZ)
    assert str(R) == "Laurent polynomial ring in x, y over ZZ with lex order"


def test_LaurentRing_new():
    x, y = symbols("x y")
    Rp, _, _ = ring("x y", ZZ)
    Rl, _, _ = laurent_ring("x y", ZZ)

    assert Rl(1) == Rl.one
    assert Rl(0) == Rl.zero
    assert Rl(Rp(1)) == Rl.one
    assert Rl(Rp(0)) == Rl.zero
    assert Rl(Rl.one) == Rl.one
    assert Rl(Rl.zero) == Rl.zero
    assert Rl(x) == Rl.x
    assert Rl(y) == Rl.y
    assert Rl((x + y)/y) == (Rl.x + Rl.y)/Rl.y

    assert Rl.new(Rp.one) == Rl.one
    assert Rl.new(Rp.x, Rp.y) == Rl.new(Rp.x) / Rl.new(Rp.y) == Rl(x/y)

    assert Rl.ground_new(1) == Rl.one
    assert Rl.ground_new(0) == Rl.zero
    assert Rl.ground_new(ZZ(1)) == Rl.one
    assert Rl.ground_new(ZZ(0)) == Rl.zero

    assert Rl.ring_new(1) == Rl.one
    assert Rl.ring_new(0) == Rl.zero
    assert Rl.ring_new(Rp(1)) == Rl.one
    assert Rl.ring_new(Rp(0)) == Rl.zero

    numer = 2*Rp.x + Rp.y
    p = Rl.ring_new(numer)
    assert p.numer == numer
    assert p.denom == Rp.one

    numer = ZZ(2)
    denom = Rp.x
    p = Rl.ring_new((numer, denom))
    p2 = Rl.new(Rp(numer), Rp(denom))
    assert p.numer == numer
    assert p.denom == denom
    assert p == p2

    p = Rl.from_expr((2*x + y)/x)
    assert p.numer == 2*Rp.x + Rp.y
    assert p.denom == Rp.x
    assert p == (2*Rl.x + Rl.y)/Rl.x

    p = {(-1, 1): 2, (0, 0): 1}
    assert Rl.from_dict(p) == 1 + 2*Rl.y/Rl.x
    p = {(0, 0): 1, (1, 1): 2}
    assert Rl.from_dict(p) == 1 + 2*Rl.x*Rl.y

    raises(ValueError, Rl.from_expr, x**y)


def test_ring_field():
    x, y = symbols("x y")
    R, _, _ = laurent_ring("x y", ZZ)
    assert R.to_domain() == ZZ.laurent_poly_ring(x, y)
    assert R.to_field() == FracField([x, y], ZZ)


def test_LaurentPolyElement_init():
    Rp, _, _ = ring("x y", ZZ)
    Rl, _, _ = laurent_ring("x y", ZZ)

    numer = 2*Rp.x + Rp.y
    denom = Rp.x
    p = LaurentPolyElement(Rl, numer, denom)
    assert p.ring == Rl
    assert p.numer == numer
    assert p.denom == denom

    Rp2, _, _ = ring("x y", QQ)

    raises(AssertionError, lambda: LaurentPolyElement(Rl, numer, 1))
    raises(AssertionError, lambda: LaurentPolyElement(Rl, Rp2.one, Rl.one))
    raises(AssertionError, lambda: LaurentPolyElement(Rl, Rl.one, Rp2.one))
    raises(AssertionError, lambda: LaurentPolyElement(Rl, Rp.one, 2*Rp.one))
    raises(AssertionError, lambda: LaurentPolyElement(Rl, Rp.one, 2*Rp.x + Rp.y))
    raises(AssertionError, lambda: LaurentPolyElement(Rl, Rp.x, Rp.x))


def test_LaurentPolyElement_eq_hash():
    Rp, _, _ = ring("x y", ZZ)
    Rl, _, _ = laurent_ring("x y", ZZ)
    p1 = (2*Rl.x + Rl.y) / Rl.x
    p2 = (2*Rl.x + Rl.y) / Rl.y
    p3 = (3*Rl.x + Rl.y) / Rl.x
    p4 = (2*Rl.x + Rl.y) / Rl.x**2
    ps = [p1, p2, p3, p4]
    for i, pi in enumerate(ps):
        for j, pj in enumerate(ps):
            if i == j:
                assert hash(pi) == hash(pj)
                assert (pi == pj) is True
                assert (pi != pj) is False
            else:
                assert (pi == pj) is False
                assert (pi != pj) is True

    p = Rl.x + Rl.y
    assert (p == p.numer) is False
    assert (p != p.numer) is True


def test_LaurentPolyElement_attributes():
    Rp, _, _ = ring("x y", ZZ)
    Rl, _, _ = laurent_ring("x y", ZZ)
    p = (2*Rl.x + Rl.y) / Rl.x
    assert p.ring == Rl
    assert p.numer == 2*Rp.x + Rp.y
    assert p.denom == Rp.x


def test_LaurentPolyElement_str():
    R, x, y = laurent_ring("x y", ZZ)
    p = (2*x + y)/x
    assert str(p) == "2 + y/x"
    assert str(R.zero) == "0"
    assert str(R.one) == "1"


def test_LaurentPolyElement_as_expr():
    x, y = symbols("x y")
    R, _, _ = laurent_ring("x y", ZZ)
    p = (2*R.x + R.y)/R.x
    assert p.as_expr() == p.as_expr_fraction() == (2*x + y)/x
    assert p.as_expr(fraction=False) == p.as_expr_add() == 2 + y/x


def test_LaurentPolyElement_properties():
    R, x, y = laurent_ring("x y", ZZ)

    assert bool(R.zero) is False
    assert bool(R.one) is True
    assert bool(x) is True

    assert R.zero.is_zero is True
    assert x.is_zero is False
    assert R.one.is_one is True
    assert x.is_one is False

    assert R.one.is_term is True
    assert R(2).is_term is True
    assert x.is_term is True
    assert (2*x).is_term is True
    assert (1/x).is_term is True
    assert (x + y).is_term is False

    assert R.zero.is_ground is True
    assert x.is_ground is False
    assert R.one.is_ground is True
    assert (2*R.one).is_ground is True

    assert R.zero.LC == ZZ.zero
    assert R.one.LC == ZZ.one
    assert x.LC == ZZ.one
    assert (2*x).LC == ZZ(2)


def test_LaurentPolyElement_arith():
    Rp, _, _ = ring("x y", ZZ)
    R, x, y = laurent_ring("x y", ZZ)
    p = (2*x + y)/x
    q = (3*x + 2*y)/y
    assert +p == p
    assert -p == (-2*x - y)/x
    assert p + q == (3*x**2 + 4*y*x + y**2) / (x*y)
    assert p - q == (-3*x**2 + y**2) / (x*y)
    assert p*q == (6*x**2 + 7*y*x + 2*y**2) / (x*y)
    assert p**2 == (4*x**2 + 4*y*x + y**2) / (x**2)
    assert p/x == (2*x + y) / x**2
    assert p/y == (2*x + y) / (x*y)
    assert (x**2 + x**2*y) / x == x + x*y
    assert (x**2 + x**2*y) / x**2 == 1 + y
    assert (x**2 + x**2*y) / x**3 == (1 + y) / x

    assert p + Rp.x == p + R.x
    assert p - Rp.x == p - R.x
    assert p * Rp.x == p * R.x
    assert p / Rp.x == p / R.x
    assert Rp.x + x == Rp.x + x
    assert Rp.x - x == Rp.x - x
    assert Rp.x * x == Rp.x * x
    assert Rp.x / x == Rp.x / x

    R2, x2, y2 = laurent_ring("x y", QQ)
    p = x
    p2 = x2
    assert p + 1 == x + R.one
    assert 1 + p == x + R.one
    assert p - 1 == x - R.one
    assert 1 - p == R.one - x
    assert p * 2 == 2*x
    assert 2 * p == 2*x
    assert p * p == x**2

    raises(TypeError, lambda: p + ())
    raises(TypeError, lambda: () + p)
    raises(TypeError, lambda: p + p2)
    raises(TypeError, lambda: p2 + p)

    raises(TypeError, lambda: p - ())
    raises(TypeError, lambda: () - p)
    raises(TypeError, lambda: p - p2)
    raises(TypeError, lambda: p2 - p)

    raises(TypeError, lambda: p * ())
    raises(TypeError, lambda: () * p)
    raises(TypeError, lambda: p * p2)
    raises(TypeError, lambda: p2 * p)

    assert x**0 == R.one
    assert x**1 == x
    assert x**2 == x*x
    assert x**-1 == 1/x
    assert x**-2 == 1/x**2
    raises(ZeroDivisionError, lambda: R.zero ** -1)
    raises(NotImplementedError, lambda: (x + y)**-1)
    raises(TypeError, lambda: p ** p)

    assert (x + 1) / x == 1 + 1/x
    assert x / x == R.one
    assert (x + y) / (x*y) == 1/x + 1/y
    raises(ZeroDivisionError, lambda: x / 0)
    raises(ZeroDivisionError, lambda: x / R.zero)
    raises(NotImplementedError, lambda: 1 / (x + y))
    raises(TypeError, lambda: x / ())
    raises(TypeError, lambda: () / x)
    raises(TypeError, lambda: p / p2)


def test_LaurentPolyElement_division():
    R, x, y = laurent_ring("x y", ZZ)
    a = (2*x + y)/x
    q = (1 + x)/y
    b = a*q
    assert b.exquo(a) == q

    raises(ExactQuotientFailed, lambda: a.exquo(b))
    raises(ZeroDivisionError, lambda: a.exquo(R.zero))

    assert ((x**2 - y**2)/x) % (x - y) == R.zero
    assert ((x**2 - y**2)/x + 1) % (x - y) == y/x
    assert R.x % 1 == R.zero
    assert 1 % R.x == R.one

    assert ((x**2 - y**2)/x) // (x - y) == (x + y)/x
    assert x // 1 == x
    assert 1 // x == R.zero

    p1 = (x**2 - y**2)/x
    p2 = (x**2 - y**2)/x + 1
    assert (p1 // p2) * p2 + (p1 % p2) == p1
    assert divmod(p1, p2) == (p1 // p2, p1 % p2)
    assert divmod(1, R.x) == divmod(R.one, R.x) == (R.zero, R.one)

    raises(ZeroDivisionError, lambda: R.one % R.zero)
    raises(ZeroDivisionError, lambda: R.one // R.zero)
    raises(ZeroDivisionError, lambda: divmod(R.one, R.zero))

    raises(TypeError, lambda: R.x % ())
    raises(TypeError, lambda: () % R.x)
    raises(TypeError, lambda: R.x // ())
    raises(TypeError, lambda: () // R.x)
    raises(TypeError, lambda: divmod(R.x, ()))
    raises(TypeError, lambda: divmod((), R.x))

    assert divmod((x**2 - y**2)/x, x - y) == ((x + y)/x, R.zero)



def test_LaurentPolyElement_gcd():
    R, x, y = laurent_ring("x y", ZZ)
    p = (2*x + y)/x
    q = (3*x + 2*y)/y
    g = (x + y)/x
    assert (p*g).gcd(q*g) == g


def test_LaurentPolyElement_factor_list():
    R, x, y = laurent_ring("x y", ZZ)
    p = (2*x + y)/x**2
    q = (3*x + 2*y)/y
    factors = [
        (2*x + y, 1),
        (3*x + 2*y, 1),
        (1/y, 1),
        (1/x, 2),
    ]
    assert (p*q)._factor_list() == (ZZ(1), factors)
