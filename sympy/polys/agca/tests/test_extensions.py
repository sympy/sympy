from sympy.polys import QQ
from sympy.polys.polytools import Poly
from sympy.polys.polyerrors import NotInvertible
from sympy.polys.agca.extensions import FiniteExtension
from sympy.testing.pytest import raises
from sympy.abc import x, t


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
    raises(NotImplementedError, lambda: A(1).inverse())

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

    # Test mod
    K = FiniteExtension(Poly(x**3 + 1))
    xf = K(x)
    assert (xf**2 - 1) % (xf - 1) == K.zero

    assert K.from_sympy(x) == xf
    assert K.to_sympy(xf) == x

    KZ = FiniteExtension(Poly(x**2 + 1, x, domain='ZZ'))
    KQ = FiniteExtension(Poly(x**2 + 1, x, domain='QQ'))
    assert KZ.set_domain(QQ) == KQ

    # Test exquo
    K = FiniteExtension(Poly(x**4 + 1))
    xf = K(x)
    assert K.exquo(xf**2 - 1, xf - 1) == xf + 1

    # Test from_MonogenicFiniteExtension
    K1 = FiniteExtension(Poly(x**2 + 1))
    K2 = QQ[x]
    x1, x2 = K1(x), K2(x)
    assert K1.convert(x2) == x1
    assert K2.convert(x1) == x2

    K = FiniteExtension(Poly(x**2 - 1, domain=QQ))
    assert K.convert_from(QQ(1, 2), QQ) == K.one/2
