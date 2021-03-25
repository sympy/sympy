from sympy.testing.pytest import raises

from sympy.polys import ZZ, QQ
from sympy.polys.matrices.domainscalar import DomainScalar


def test_DomainScalar_init():
    A = DomainScalar(ZZ(1), ZZ)
    assert A.element ==  ZZ(1)
    assert A.domain == ZZ


def test_DomainScalar_new():
    A = DomainScalar(ZZ(1), ZZ)
    B = A.new(ZZ(4), ZZ)
    assert B == DomainScalar(ZZ(4), ZZ)


def test_DomainScalar_repr():
    A = DomainScalar(ZZ(1), ZZ)
    assert A.__repr__() == '1'


def test_DomainScalar_from_sympy():
    expr = ZZ(1)
    B = DomainScalar.from_sympy(expr)
    assert B == DomainScalar(ZZ(1), ZZ)


def test_DomainScalar_to_domain():
    A = DomainScalar(ZZ(1), ZZ)
    B = A.to_domain(QQ)
    assert B == DomainScalar(QQ(1), QQ)


def test_DomainScalar_unify():
    A = DomainScalar(ZZ(1), ZZ)
    B = DomainScalar(QQ(2), QQ)
    A, B = A.unify(B)
    assert A.domain == B.domain == QQ


def test_DomainScalar_add():
    A = DomainScalar(ZZ(1), ZZ)
    B = DomainScalar(QQ(2), QQ)
    C = A.__add__(B)
    assert C == DomainScalar(QQ(3), QQ)


def test_DomainScalar_sub():
    A = DomainScalar(ZZ(1), ZZ)
    B = DomainScalar(QQ(2), QQ)
    C = A.__sub__(B)
    assert C == DomainScalar(QQ(-1), QQ)


def test_DomainScalar_mul():
    A = DomainScalar(ZZ(1), ZZ)
    B = DomainScalar(QQ(2), QQ)
    C = A.__mul__(B)
    assert C == DomainScalar(QQ(2), QQ)


def test_DomainScalar_floordiv():
    A = DomainScalar(ZZ(-5), ZZ)
    B = DomainScalar(QQ(2), QQ)
    C = A.__floordiv__(B)
    assert C == DomainScalar(QQ(-5)/QQ(2), QQ)


def test_DomainScalar_pow():
    A = DomainScalar(ZZ(-5), ZZ)
    B = A.__pow__(2)
    assert B == DomainScalar(ZZ(25), ZZ)

    raises(ZeroDivisionError, lambda: DomainScalar(ZZ(0), ZZ).__pow__(-2))


def test_DomainScalar_pos():
    A = DomainScalar(ZZ(-5), ZZ)
    B = A.__pos__()
    assert B == DomainScalar(ZZ(5), ZZ)
