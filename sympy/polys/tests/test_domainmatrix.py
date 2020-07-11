from sympy.core.compatibility import HAS_GMPY
from sympy.polys import ZZ, QQ
from sympy.polys.domainmatrix import (
        DDM,
        DDMBadInputError, DDMDomainError, DDMShapeError,
        ddm_add, ddm_sub, ddm_neg, ddm_mmul,
        )

from sympy.testing.pytest import raises


def test_DDM_init():
    items = [[ZZ(0), ZZ(1), ZZ(2)], [ZZ(3), ZZ(4), ZZ(5)]]
    shape = (2, 3)
    ddm = DDM(items, shape, ZZ)
    assert ddm.shape == shape
    assert ddm.rows == 2
    assert ddm.cols == 3
    assert ddm.domain == ZZ

    raises(DDMBadInputError, lambda: DDM([[ZZ(2), ZZ(3)]], (2, 2), ZZ))
    raises(DDMBadInputError, lambda: DDM([[ZZ(1)], [ZZ(2), ZZ(3)]], (2, 2), ZZ))


def test_DDM_str():
    ddm = DDM([[ZZ(0), ZZ(1)], [ZZ(2), ZZ(3)]], (2, 2), ZZ)
    if HAS_GMPY:
        assert str(ddm) == 'DDM([[mpz(0), mpz(1)], [mpz(2), mpz(3)]], (2, 2), ZZ)'
    else:
        assert str(ddm) == 'DDM([[0, 1], [2, 3]], (2, 2), ZZ)'


def test_DDM_eq():
    items = [[ZZ(0), ZZ(1)], [ZZ(2), ZZ(3)]]
    ddm1 = DDM(items, (2, 2), ZZ)
    ddm2 = DDM(items, (2, 2), ZZ)

    assert (ddm1 == ddm1) is True
    assert (ddm1 == items) is False
    assert (items == ddm1) is False
    assert (ddm1 == ddm2) is True
    assert (ddm2 == ddm1) is True

    assert (ddm1 != ddm1) is False
    assert (ddm1 != items) is True
    assert (items != ddm1) is True
    assert (ddm1 != ddm2) is False
    assert (ddm2 != ddm1) is False

    ddm3 = DDM([[ZZ(0), ZZ(1)], [ZZ(3), ZZ(3)]], (2, 2), ZZ)
    ddm3 = DDM(items, (2, 2), QQ)

    assert (ddm1 == ddm3) is False
    assert (ddm3 == ddm1) is False
    assert (ddm1 != ddm3) is True
    assert (ddm3 != ddm1) is True


def test_DDM_zeros():
    ddmz = DDM.zeros((3, 4), QQ)
    assert list(ddmz) == [[QQ(0)] * 4] * 3
    assert ddmz.shape == (3, 4)
    assert ddmz.domain == QQ


def test_DDM_copy():
    ddm1 = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    ddm2 = ddm1.copy()
    assert (ddm1 == ddm2) is True
    ddm1[0][0] = QQ(-1)
    assert (ddm1 == ddm2) is False
    ddm2[0][0] = QQ(-1)
    assert (ddm1 == ddm2) is True


def test_DDM_add():
    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    B = DDM([[ZZ(3)], [ZZ(4)]], (2, 1), ZZ)
    C = DDM([[ZZ(4)], [ZZ(6)]], (2, 1), ZZ)
    AQ = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    assert A + B == A.add(B) == C
    raises(DDMShapeError, lambda: A + DDM([[ZZ(5)]], (1, 1), ZZ))
    raises(TypeError, lambda: A + ZZ(1))
    raises(TypeError, lambda: ZZ(1) + A)
    raises(DDMDomainError, lambda: A + AQ)
    raises(DDMDomainError, lambda: AQ + A)


def test_DDM_sub():
    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    B = DDM([[ZZ(3)], [ZZ(4)]], (2, 1), ZZ)
    C = DDM([[ZZ(-2)], [ZZ(-2)]], (2, 1), ZZ)
    AQ = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    D = DDM([[ZZ(5)]], (1, 1), ZZ)
    assert A - B == A.sub(B) == C
    raises(TypeError, lambda: A - ZZ(1))
    raises(TypeError, lambda: ZZ(1) - A)
    raises(DDMShapeError, lambda: A - D)
    raises(DDMShapeError, lambda: D - A)
    raises(DDMShapeError, lambda: A.sub(D))
    raises(DDMShapeError, lambda: D.sub(A))
    raises(DDMDomainError, lambda: A - AQ)
    raises(DDMDomainError, lambda: AQ - A)
    raises(DDMDomainError, lambda: A.sub(AQ))
    raises(DDMDomainError, lambda: AQ.sub(A))


def test_DDM_neg():
    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    An = DDM([[ZZ(-1)], [ZZ(-2)]], (2, 1), ZZ)
    assert -A == A.neg() == An
    assert -An == An.neg() == A


def test_DDM_mul():
    A = DDM([[ZZ(1)]], (1, 1), ZZ)
    raises(TypeError, lambda: [[1]] * A)
    raises(TypeError, lambda: A * [[1]])


def test_DDM_mmul():
    # ZZZ: Use a @ b?
    # ZZZ: test 0xn matrices
    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    B = DDM([[ZZ(3), ZZ(4)]], (1, 2), ZZ)
    AB = DDM([[ZZ(3), ZZ(4)], [ZZ(6), ZZ(8)]], (2, 2), ZZ)
    BA = DDM([[ZZ(11)]], (1, 1), ZZ)
    Bq = DDM([[QQ(3), QQ(4)]], (1, 2), QQ)
    C = DDM([[ZZ(1)]], (1, 1), ZZ)
    assert A * B == A.mul(B) == AB
    assert B * A == B.mul(A) == BA
    raises(DDMDomainError, lambda: A * Bq)
    raises(DDMDomainError, lambda: Bq * A)
    assert A * C == A.mul(C) == A
    raises(DDMShapeError, lambda: C * A)
    raises(DDMShapeError, lambda: C.mul(A))


def test_ddm_add():
    a = [[1, 2], [3, 4]]
    b = [[5, 6], [7, 8]]
    ddm_add(a, b)
    assert a == [[6, 8], [10, 12]]


def test_ddm_sub():
    a = [[1, 2], [3, 4]]
    b = [[5, 6], [7, 8]]
    ddm_sub(a, b)
    assert a == [[-4, -4], [-4, -4]]

def test_ddm_neg():
    a = [[1, 2], [3, 4]]
    ddm_neg(a)
    assert a == [[-1, -2], [-3, -4]]


def test_ddm_mmul():
    a = [[1, 2, 3], [4, 5, 6]]
    b = [[1, 2], [3, 4], [5, 6]]
    c1 = [[0, 0], [0, 0]]
    c2 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    ddm_mmul(a, b, c1)
    assert c1 == [[22, 28], [49, 64]]
    ddm_mmul(b, a, c2)
    assert c2 == [[9, 12, 15], [19, 26, 33], [29, 40, 51]]
    b3 = [[1], [2], [3]]
    c3 = [[0], [0]]
    ddm_mmul(a, b3, c3)
    assert c3 == [[14], [32]]
