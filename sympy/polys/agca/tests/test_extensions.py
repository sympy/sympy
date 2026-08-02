from __future__ import annotations
from sympy.core.symbol import symbols
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.polys import QQ, ZZ
from sympy.polys.domains.ring import Ring
from sympy.polys.domains.ringextension import RingExtension
from sympy.polys.polytools import Poly
from sympy.polys.polyerrors import (CoercionFailed, DomainError, NotInvertible,
    NotReversible)
from sympy.polys.agca.extensions import ExtensionElement, FiniteExtension
from sympy.polys.domainmatrix import DomainMatrix

from sympy.testing.pytest import raises

from sympy.abc import x, y, t


def test_FiniteExtension():
    # Gaussian integers
    A = FiniteExtension(Poly(x**2 + 1, x))
    assert A.rank == 2
    assert str(A) == 'ZZ[x]/(x**2 + 1)'
    i = A.generator
    assert i.parent() is A

    assert i*i == A(-1)
    raises(TypeError, lambda: i*())

    assert A.basis == (A.one, i)
    assert A(1) == A.one
    assert i**2 == A(-1)
    assert i**2 != -1  # no coercion
    assert (2 + i)*(1 - i) == 3 - i
    assert (1 + i)**8 == A(16)
    assert A(1).inverse() == A(1)
    raises(NotImplementedError, lambda: A(2).inverse())

    # Finite field of order 27
    F = FiniteExtension(Poly(x**3 - x + 1, x, modulus=3))
    assert F.rank == 3
    a = F.generator  # also generates the cyclic group F - {0}
    assert F.basis == (F(1), a, a**2)
    assert a**27 == a
    assert a**26 == F(1)
    assert a**13 == F(-1)
    assert a**9 == a + 1
    assert a**3 == a - 1
    assert a**6 == a**2 + a + 1
    assert F(x**2 + x).inverse() == 1 - a
    assert F(x + 2)**(-1) == F(x + 2).inverse()
    assert a**19 * a**(-19) == F(1)
    assert (a - 1) / (2*a**2 - 1) == a**2 + 1
    assert (a - 1) // (2*a**2 - 1) == a**2 + 1
    assert 2/(a**2 + 1) == a**2 - a + 1
    assert (a**2 + 1)/2 == -a**2 - 1
    raises(NotInvertible, lambda: F(0).inverse())

    # Function field of an elliptic curve
    K = FiniteExtension(Poly(t**2 - x**3 - x + 1, t, field=True))
    assert K.rank == 2
    assert str(K) == 'ZZ(x)[t]/(t**2 - x**3 - x + 1)'
    y = K.generator
    c = 1/(x**3 - x**2 + x - 1)
    assert ((y + x)*(y - x)).inverse() == K(c)
    assert (y + x)*(y - x)*c == K(1)  # explicit inverse of y + x


def test_FiniteExtension_eq_hash():
    # Test eq and hash
    p1 = Poly(x**2 - 2, x, domain=ZZ)
    p2 = Poly(x**2 - 2, x, domain=QQ)
    K1 = FiniteExtension(p1)
    K2 = FiniteExtension(p2)
    assert K1 == FiniteExtension(Poly(x**2 - 2))
    assert K2 != FiniteExtension(Poly(x**2 - 2))
    assert len({K1, K2, FiniteExtension(p1)}) == 2


def test_FiniteExtension_mod():
    # Test mod
    K = FiniteExtension(Poly(x**3 + 1, x, domain=QQ))
    xf = K(x)
    assert (xf**2 - 1) % 1 == K.zero
    raises(ZeroDivisionError, lambda: 1 % (xf**2 - 1))
    assert (xf**2 - 1) / (xf - 1) == xf + 1
    assert (xf**2 - 1) // (xf - 1) == xf + 1
    assert (xf**2 - 1) % (xf - 1) == K.zero
    raises(ZeroDivisionError, lambda: (xf**2 - 1) % 0)
    raises(TypeError, lambda: xf % [])
    raises(TypeError, lambda: [] % xf)

    # Test mod over ring
    K = FiniteExtension(Poly(x**3 + 1, x, domain=ZZ))
    xf = K(x)
    assert (xf**2 - 1) % 1 == K.zero
    raises(NotImplementedError, lambda: (xf**2 - 1) % (xf - 1))


def test_FiniteExtension_from_sympy():
    # Test to_sympy/from_sympy
    K = FiniteExtension(Poly(x**3 + 1, x, domain=ZZ))
    xf = K(x)
    assert K.from_sympy(x) == xf
    assert K.to_sympy(xf) == x


def test_FiniteExtension_set_domain():
    KZ = FiniteExtension(Poly(x**2 + 1, x, domain='ZZ'))
    KQ = FiniteExtension(Poly(x**2 + 1, x, domain='QQ'))
    assert KZ.set_domain(QQ) == KQ


def test_FiniteExtension_ring_extension():
    K = FiniteExtension(Poly(x**3 + 2*x + 1, x, domain=QQ))
    xK = K.generator

    assert isinstance(K, Ring)
    assert isinstance(K, RingExtension)
    assert K.dom == QQ
    assert K.gens == (xK,)
    assert K.ngens == 1
    assert K.to_dict(3*xK**2 + 2) == {(2,): QQ(3), (0,): QQ(2)}


def test_FiniteExtension_modulus_degree():
    raises(ValueError, lambda: FiniteExtension(Poly(0, x, domain=QQ)))
    raises(ValueError, lambda: FiniteExtension(Poly(1, x, domain=QQ)))


def test_FiniteExtension_field_contract():
    K = FiniteExtension(Poly(x**2 - 2, x, domain=QQ))
    xK = K.generator

    assert K.modulus_is_irreducible is True
    assert K.is_Field is True
    assert K.is_Ring is True
    assert K.is_PID is True
    assert K.has_assoc_Field is True
    assert K.has_assoc_Ring is False
    assert K.get_field() is K
    raises(DomainError, K.get_ring)

    assert K.gcd(K.zero, K.zero) == K.one
    assert K.lcm(xK, xK + 1) == xK*(xK + 1)
    s, t, h = K.gcdex(xK, xK + 1)
    assert s*xK + t*(xK + 1) == h == K.one
    assert K.div(xK + 1, xK) == ((xK + 1)/xK, K.zero)
    assert K.rem(xK + 1, xK) == K.zero
    assert K.revert(xK) == xK/2
    assert K.numer(xK) == xK
    assert K.denom(xK) == K.one
    assert K.from_sympy(1/x) == xK/2
    assert K.to_sympy(K.from_sympy((x + 1)/(x - 1))) == 2*x + 3

    matrix = DomainMatrix([[xK, K.one]], (1, 2), K)
    assert matrix.nullspace() == DomainMatrix([[-K.one, xK]], (1, 2), K)

    F = FiniteExtension(Poly(x**3 - x + 1, x, modulus=3))
    assert F.modulus_is_irreducible is True
    assert F.is_Field is True

    R = FiniteExtension(Poly(x**2 - 1, x, domain=QQ))
    xR = R.generator
    assert R.modulus_is_irreducible is False
    assert R.is_Field is False
    assert R.is_PID is False
    assert R.has_assoc_Field is False
    assert R.has_assoc_Ring is True
    assert R.get_ring() is R
    raises(DomainError, R.get_field)
    raises(NotImplementedError, lambda: R.gcd(xR, R.one))
    assert R.numer(xR) == xR
    assert R.denom(xR) == R.one
    assert R.is_unit(xR) is True
    assert R.revert(xR) == xR
    assert R.one % xR == R.zero
    assert R.is_unit(xR - 1) is False
    zero_divisor = xR - 1
    raises(ZeroDivisionError, lambda: R.one % zero_divisor)
    raises(ZeroDivisionError, lambda: R.rem(R.one, zero_divisor))
    assert zero_divisor % zero_divisor == R.zero
    assert R.div(zero_divisor, zero_divisor) == (
        R.quo(zero_divisor, zero_divisor),
        R.rem(zero_divisor, zero_divisor),
    ) == (R.one, R.zero)
    assert R.zero % zero_divisor == R.zero
    assert R.div(R.zero, zero_divisor) == (
        R.quo(R.zero, zero_divisor),
        R.rem(R.zero, zero_divisor),
    ) == (R.zero, R.zero)
    raises(NotReversible, lambda: R.revert(zero_divisor))

    S = FiniteExtension(Poly(y**2 - 1, y, domain=R))
    assert S.modulus_is_irreducible is None
    assert S.is_Field is False


def test_FiniteExtension_unknown_irreducibility():
    from sympy.polys.domains import GF

    K = GF(5).frac_field(x)
    A = FiniteExtension(Poly(y**2 - x, y, domain=K))

    assert A.modulus_is_irreducible is None
    assert A.is_Field is False
    assert A.has_assoc_Field is False
    assert A.has_assoc_Ring is True
    assert A.generator**2 == A(x)
    raises(NotImplementedError, lambda: A.gcd(A.one, A.one))

    B = FiniteExtension(Poly(y**2 - x**2, y, domain=K))
    assert B.modulus_is_irreducible is None
    assert B.is_Field is False
    assert (B.generator - B(x))*(B.generator + B(x)) == B.zero


def test_FiniteExtension_irreducibility_checked_once():
    class CountingPoly(Poly):
        irreducibility_checks = 0

        @property
        def is_irreducible(self):
            type(self).irreducibility_checks += 1
            return super().is_irreducible

    K = FiniteExtension(CountingPoly(x**2 - 2, x, domain=QQ))

    assert K.modulus_is_irreducible is True
    assert K.is_Field is True
    assert K.modulus_is_irreducible is True
    assert CountingPoly.irreducibility_checks == 1


def test_FiniteExtension_exquo():
    # Test exquo
    K = FiniteExtension(Poly(x**4 + 1))
    xf = K(x)
    assert K.exquo(xf**2 - 1, xf - 1) == xf + 1

    K1 = FiniteExtension(Poly(x**3 - 101*x**2 - 1141*x - 15066, x, domain=QQ))
    x1 = K1(x)
    assert K1.exquo(-25*x1**2 - 1797*x1 - 15106, 10 - x1) == x1**2 - 66*x1 - 4


def test_FiniteExtension_convert():
    # Test from_MonogenicFiniteExtension
    K1 = FiniteExtension(Poly(x**2 + 1))
    K2 = QQ[x]
    x1, x2 = K1(x), K2(x)
    assert K1.convert(x2) == x1
    assert K2.convert(x1) == x2

    K = FiniteExtension(Poly(x**2 - 1, domain=QQ))
    assert K.convert_from(QQ(1, 2), QQ) == K.one/2


def test_FiniteExtension_convert_provenance():
    K = FiniteExtension(Poly(x**2 - 2, x, domain=QQ))
    malformed = ExtensionElement(K.ring.convert(x**2), K)

    assert malformed != K(2)
    assert K.convert(malformed) == K(2)
    assert K.convert_from(malformed, K) == K(2)

    K2 = FiniteExtension(Poly(t**2 - 3, t, domain=K))
    assert K2.convert(malformed) == K2(2)
    assert K2.convert_from(malformed, K) == K2(2)

    other = FiniteExtension(Poly(x**2 - 3, x, domain=QQ))
    raises(CoercionFailed, lambda: K.convert(other.generator))
    raises(CoercionFailed, lambda: K.convert_from(other.generator, other))
    raises(CoercionFailed, lambda: K2.convert_from(other.generator, other))


def test_FiniteExtension_division_ring():
    # Test division in FiniteExtension over a ring
    KQ = FiniteExtension(Poly(x**2 - 1, x, domain=QQ))
    KZ = FiniteExtension(Poly(x**2 - 1, x, domain=ZZ))
    KQt = FiniteExtension(Poly(x**2 - 1, x, domain=QQ[t]))
    KQtf = FiniteExtension(Poly(x**2 - 1, x, domain=QQ.frac_field(t)))
    assert KQ.is_Field is False
    assert KZ.is_Field is False
    assert KQt.is_Field is False
    assert KQtf.is_Field is False
    for K in KQ, KZ, KQt, KQtf:
        xK = K.convert(x)
        assert xK / K.one == xK
        assert xK // K.one == xK
        assert xK % K.one == K.zero
        raises(ZeroDivisionError, lambda: xK / K.zero)
        raises(ZeroDivisionError, lambda: xK // K.zero)
        raises(ZeroDivisionError, lambda: xK % K.zero)
        if K.domain.is_Field:
            assert xK / xK == K.one
            assert xK // xK == K.one
            assert xK % xK == K.zero
        else:
            raises(NotImplementedError, lambda: xK / xK)
            raises(NotImplementedError, lambda: xK // xK)
            raises(NotImplementedError, lambda: xK % xK)


def test_FiniteExtension_composite_domain_embedding():
    u = symbols('u')
    K = QQ.frac_field(x)
    L = FiniteExtension(Poly(y**2 - x, y, domain=K))
    yL = L.generator
    R = L.poly_ring(u)
    F = L.frac_field(u)

    assert L.is_Field is True

    yR = R.convert(yL)
    yF = F.convert(yL)
    assert yR == R.ring.ground_new(yL)
    assert yF == F.field.ground_new(yL)

    uR = R.gens[0]
    uF = F.gens[0]
    assert (uR + yR)*(uR - yR) == uR**2 - R.from_sympy(x)

    poly = uR + yR
    assert R.convert(F.convert(poly, R), F) == poly

    value = F.from_sympy((u + y)/(u - y))
    assert value == (uF + yF)/(uF - yF)
    assert (uF - yF)*value == uF + yF
    assert F.to_sympy(value) == (u + y)/(u - y)

    matrix = DomainMatrix([[uF + yF, F.one]], (1, 2), F)
    expected = DomainMatrix([[-F.one, uF + yF]], (1, 2), F)
    assert matrix.nullspace() == expected

    other = FiniteExtension(Poly(y**2 - x - 1, y, domain=K))
    raises(CoercionFailed, lambda: R.convert(other.generator))
    raises(CoercionFailed, lambda: F.convert(other.generator))


def test_FiniteExtension_Poly():
    K = FiniteExtension(Poly(x**2 - 2))
    p = Poly(x, y, domain=K)
    assert p.domain == K
    assert p.as_expr() == x
    assert (p**2).as_expr() == 2

    K = FiniteExtension(Poly(x**2 - 2, x, domain=QQ))
    # t**2 - 3 is irreducible over QQ(sqrt(2)), but factoring over a
    # MonogenicFiniteExtension is not implemented yet.
    K2 = FiniteExtension(Poly(t**2 - 3, t, domain=K))
    assert K2.modulus_is_irreducible is None
    assert K2.is_Field is False
    assert str(K2) == 'QQ[x]/(x**2 - 2)[t]/(t**2 - 3)'
    assert K2.convert(K.generator) == K2.convert(x)

    other = FiniteExtension(Poly(y**2 - 3, y, domain=QQ))
    raises(CoercionFailed, lambda: K2.convert(other.generator))

    eK = K2.convert(x + t)
    assert K2.to_sympy(eK) == x + t
    assert K2.to_sympy(eK ** 2) == 5 + 2*x*t
    p = Poly(x + t, y, domain=K2)
    assert p**2 == Poly(5 + 2*x*t, y, domain=K2)


def test_FiniteExtension_sincos_jacobian():
    # Use FiniteExtensino to compute the Jacobian of a matrix involving sin
    # and cos of different symbols.
    r, p, t = symbols('rho, phi, theta')
    elements = [
        [sin(p)*cos(t), r*cos(p)*cos(t), -r*sin(p)*sin(t)],
        [sin(p)*sin(t), r*cos(p)*sin(t),  r*sin(p)*cos(t)],
        [       cos(p),       -r*sin(p),                0],
    ]

    def make_extension(K):
        K = FiniteExtension(Poly(sin(p)**2+cos(p)**2-1, sin(p), domain=K[cos(p)]))
        K = FiniteExtension(Poly(sin(t)**2+cos(t)**2-1, sin(t), domain=K[cos(t)]))
        return K

    Ksc1 = make_extension(ZZ[r])
    Ksc2 = make_extension(ZZ)[r]

    for K in [Ksc1, Ksc2]:
        elements_K = [[K.convert(e) for e in row] for row in elements]
        J = DomainMatrix(elements_K, (3, 3), K)
        det = J.charpoly()[-1] * (-K.one)**3
        assert det == K.convert(r**2*sin(p))

def test_nullspace_finiteextension():
    mat = DomainMatrix({0: {0: 10, 1: 12, 2: 44}, 1: {0: 47, 1: 56, 2: 67},
                        2: {0: 22, 1: 37, 2: 35}}, (3, 3), QQ)
    field = FiniteExtension(Poly(x**3 - 101*x**2 - 1141*x - 15066, x, domain=QQ))
    expected_result = DomainMatrix({0: {i: field.convert(e)
            for i, e in enumerate([44*x - 1660, 67*x + 1398, x**2 - 66*x - 4])}}, (1, 3), field)
    result = (mat.convert_to(field) - field(x) * mat.eye(mat.shape, field)).nullspace()
    assert result == expected_result
