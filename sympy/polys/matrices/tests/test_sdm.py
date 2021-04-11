"""
Tests for the basic functionality of the SDM class.
"""

from sympy.core.compatibility import HAS_GMPY
from sympy.testing.pytest import raises

from sympy import QQ, ZZ
from sympy.polys.matrices.sdm import SDM
from sympy.polys.matrices.ddm import DDM
from sympy.polys.matrices.exceptions import (DDMBadInputError, DDMDomainError,
        DDMShapeError)


def test_SDM():
    A = SDM({0:{0:ZZ(1)}}, (2, 2), ZZ)
    assert A.domain == ZZ
    assert A.shape == (2, 2)
    assert dict(A) == {0:{0:ZZ(1)}}

    raises(DDMBadInputError, lambda: SDM({5:{1:ZZ(0)}}, (2, 2), ZZ))
    raises(DDMBadInputError, lambda: SDM({0:{5:ZZ(0)}}, (2, 2), ZZ))


def test_DDM_str():
    sdm = SDM({0:{0:ZZ(1)}, 1:{1:ZZ(1)}}, (2, 2), ZZ)
    assert str(sdm) == '{0: {0: 1}, 1: {1: 1}}'
    if HAS_GMPY: # pragma: no cover
        assert repr(sdm) == 'SDM({0: {0: mpz(1)}, 1: {1: mpz(1)}}, (2, 2), ZZ)'
    else:        # pragma: no cover
        assert repr(sdm) == 'SDM({0: {0: 1}, 1: {1: 1}}, (2, 2), ZZ)'


def test_SDM_new():
    A = SDM({0:{0:ZZ(1)}}, (2, 2), ZZ)
    B = A.new({}, (2, 2), ZZ)
    assert B == SDM({}, (2, 2), ZZ)


def test_SDM_copy():
    A = SDM({0:{0:ZZ(1)}}, (2, 2), ZZ)
    B = A.copy()
    assert A == B
    A[0][0] = ZZ(2)
    assert A != B


def test_SDM_from_list():
    A = SDM.from_list([[ZZ(0), ZZ(1)], [ZZ(1), ZZ(0)]], (2, 2), ZZ)
    assert A == SDM({0:{1:ZZ(1)}, 1:{0:ZZ(1)}}, (2, 2), ZZ)

    raises(DDMBadInputError, lambda: SDM.from_list([[ZZ(0)], [ZZ(0), ZZ(1)]], (2, 2), ZZ))
    raises(DDMBadInputError, lambda: SDM.from_list([[ZZ(0), ZZ(1)]], (2, 2), ZZ))


def test_SDM_to_list():
    A = SDM({0:{1: ZZ(1)}}, (2, 2), ZZ)
    assert A.to_list() == [[ZZ(0), ZZ(1)], [ZZ(0), ZZ(0)]]

    A = SDM({}, (0, 2), ZZ)
    assert A.to_list() == []

    A = SDM({}, (2, 0), ZZ)
    assert A.to_list() == [[], []]


def test_SDM_from_ddm():
    A = DDM([[ZZ(1), ZZ(0)], [ZZ(1), ZZ(0)]], (2, 2), ZZ)
    B = SDM.from_ddm(A)
    assert B.domain == ZZ
    assert B.shape == (2, 2)
    assert dict(B) == {0:{0:ZZ(1)}, 1:{0:ZZ(1)}}


def test_SDM_to_ddm():
    A = SDM({0:{1: ZZ(1)}}, (2, 2), ZZ)
    B = DDM([[ZZ(0), ZZ(1)], [ZZ(0), ZZ(0)]], (2, 2), ZZ)
    assert A.to_ddm() == B


def test_SDM_zeros():
    A = SDM.zeros((2, 2), ZZ)
    assert A.domain == ZZ
    assert A.shape == (2, 2)
    assert dict(A) == {}


def test_SDM_eye():
    A = SDM.eye(2, ZZ)
    assert A.domain == ZZ
    assert A.shape == (2, 2)
    assert dict(A) == {0:{0:ZZ(1)}, 1:{1:ZZ(1)}}


def test_SDM_transpose():
    A = SDM({0:{0:ZZ(1), 1:ZZ(2)}, 1:{0:ZZ(3), 1:ZZ(4)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(1), 1:ZZ(3)}, 1:{0:ZZ(2), 1:ZZ(4)}}, (2, 2), ZZ)
    assert A.transpose() == B

    A = SDM({0:{1:ZZ(2)}}, (2, 2), ZZ)
    B = SDM({1:{0:ZZ(2)}}, (2, 2), ZZ)
    assert A.transpose() == B

    A = SDM({0:{1:ZZ(2)}}, (1, 2), ZZ)
    B = SDM({1:{0:ZZ(2)}}, (2, 1), ZZ)
    assert A.transpose() == B


def test_SDM_mul():
    A = SDM({0:{0:ZZ(2)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(4)}}, (2, 2), ZZ)
    assert A*ZZ(2) == B
    assert ZZ(2)*A == B

    raises(TypeError, lambda: A*QQ(1, 2))
    raises(TypeError, lambda: QQ(1, 2)*A)


def test_SDM_matmul():
    A = SDM({0:{0:ZZ(2)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(4)}}, (2, 2), ZZ)
    assert A.matmul(A) == B

    C = SDM({0:{0:ZZ(2)}}, (2, 2), QQ)
    raises(DDMDomainError, lambda: A.matmul(C))

    A = SDM({0:{0:ZZ(1), 1:ZZ(2)}, 1:{0:ZZ(3), 1:ZZ(4)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(7), 1:ZZ(10)}, 1:{0:ZZ(15), 1:ZZ(22)}}, (2, 2), ZZ)
    assert A.matmul(A) == B

    A22 = SDM({0:{0:ZZ(4)}}, (2, 2), ZZ)
    A32 = SDM({0:{0:ZZ(2)}}, (3, 2), ZZ)
    A23 = SDM({0:{0:ZZ(4)}}, (2, 3), ZZ)
    A33 = SDM({0:{0:ZZ(8)}}, (3, 3), ZZ)
    A22 = SDM({0:{0:ZZ(8)}}, (2, 2), ZZ)
    assert A32.matmul(A23) == A33
    assert A23.matmul(A32) == A22
    # XXX: @ not supported by SDM...
    #assert A32.matmul(A23) == A32 @ A23 == A33
    #assert A23.matmul(A32) == A23 @ A32 == A22
    #raises(DDMShapeError, lambda: A23 @ A22)
    raises(DDMShapeError, lambda: A23.matmul(A22))


def test_SDM_add():
    A = SDM({0:{1:ZZ(1)}, 1:{0:ZZ(2), 1:ZZ(3)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(1)}, 1:{0:ZZ(-2), 1:ZZ(3)}}, (2, 2), ZZ)
    C = SDM({0:{0:ZZ(1), 1:ZZ(1)}, 1:{1:ZZ(6)}}, (2, 2), ZZ)
    assert A.add(B) == C

    A = SDM({0:{1:ZZ(1)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(1)}, 1:{0:ZZ(-2), 1:ZZ(3)}}, (2, 2), ZZ)
    C = SDM({0:{0:ZZ(1), 1:ZZ(1)}, 1:{0:ZZ(-2), 1:ZZ(3)}}, (2, 2), ZZ)
    assert A.add(B) == C
    assert B.add(A) == C


def test_SDM_sub():
    A = SDM({0:{1:ZZ(1)}, 1:{0:ZZ(2), 1:ZZ(3)}}, (2, 2), ZZ)
    B = SDM({0:{0:ZZ(1)}, 1:{0:ZZ(-2), 1:ZZ(3)}}, (2, 2), ZZ)
    C = SDM({0:{0:ZZ(-1), 1:ZZ(1)}, 1:{0:ZZ(4)}}, (2, 2), ZZ)
    assert A.sub(B) == C


def test_SDM_neg():
    A = SDM({0:{1:ZZ(1)}, 1:{0:ZZ(2), 1:ZZ(3)}}, (2, 2), ZZ)
    B = SDM({0:{1:ZZ(-1)}, 1:{0:ZZ(-2), 1:ZZ(-3)}}, (2, 2), ZZ)
    assert A.neg() == B


def test_SDM_convert_to():
    A = SDM({0:{1:ZZ(1)}, 1:{0:ZZ(2), 1:ZZ(3)}}, (2, 2), ZZ)
    B = SDM({0:{1:QQ(1)}, 1:{0:QQ(2), 1:QQ(3)}}, (2, 2), QQ)
    C = A.convert_to(QQ)
    assert C == B
    assert C.domain == QQ

    D = A.convert_to(ZZ)
    assert D == A
    assert D.domain == ZZ


def test_SDM_hstack():
    A = SDM({0:{1:ZZ(1)}}, (2, 2), ZZ)
    B = SDM({1:{1:ZZ(1)}}, (2, 2), ZZ)
    AA = SDM({0:{1:ZZ(1), 3:ZZ(1)}}, (2, 4), ZZ)
    AB = SDM({0:{1:ZZ(1)}, 1:{3:ZZ(1)}}, (2, 4), ZZ)
    assert SDM.hstack(A) == A
    assert SDM.hstack(A, A) == AA
    assert SDM.hstack(A, B) == AB


def test_SDM_inv():
    A = SDM({0:{0:QQ(1), 1:QQ(2)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    B = SDM({0:{0:QQ(-2), 1:QQ(1)}, 1:{0:QQ(3, 2), 1:QQ(-1, 2)}}, (2, 2), QQ)
    assert A.inv() == B


def test_SDM_det():
    A = SDM({0:{0:QQ(1), 1:QQ(2)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    assert A.det() == QQ(-2)


def test_SDM_lu():
    A = SDM({0:{0:QQ(1), 1:QQ(2)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    L = SDM({0:{0:QQ(1)}, 1:{0:QQ(3), 1:QQ(1)}}, (2, 2), QQ)
    #U = SDM({0:{0:QQ(1), 1:QQ(2)}, 1:{0:QQ(3), 1:QQ(-2)}}, (2, 2), QQ)
    #swaps = []
    # This doesn't quite work. U has some nonzero elements in the lower part.
    #assert A.lu() == (L, U, swaps)
    assert A.lu()[0] == L


def test_SDM_lu_solve():
    A = SDM({0:{0:QQ(1), 1:QQ(2)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    b = SDM({0:{0:QQ(1)}, 1:{0:QQ(2)}}, (2, 1), QQ)
    x = SDM({1:{0:QQ(1, 2)}}, (2, 1), QQ)
    assert A.matmul(x) == b
    assert A.lu_solve(b) == x


def test_SDM_charpoly():
    A = SDM({0:{0:ZZ(1), 1:ZZ(2)}, 1:{0:ZZ(3), 1:ZZ(4)}}, (2, 2), ZZ)
    assert A.charpoly() == [ZZ(1), ZZ(-5), ZZ(-2)]


def test_SDM_nullspace():
    A = SDM({0:{0:QQ(1), 1:QQ(1)}}, (2, 2), QQ)
    assert A.nullspace()[0] == SDM({0:{0:QQ(-1), 1:QQ(1)}}, (1, 2), QQ)


def test_SDM_rref():
    eye2 = SDM({0:{0:QQ(1)}, 1:{1:QQ(1)}}, (2, 2), QQ)

    A = SDM({0:{0:QQ(1), 1:QQ(2)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    assert A.rref() == (eye2, [0, 1])

    A = SDM({0:{0:QQ(1)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    assert A.rref() == (eye2, [0, 1])

    A = SDM({0:{1:QQ(2)}, 1:{0:QQ(3), 1:QQ(4)}}, (2, 2), QQ)
    assert A.rref() == (eye2, [0, 1])

    A = SDM({0:{0:QQ(1), 1:QQ(2), 2:QQ(3)},
             1:{0:QQ(4), 1:QQ(5), 2:QQ(6)},
             2:{0:QQ(7), 1:QQ(8), 2:QQ(9)} }, (3, 3), QQ)
    Arref = SDM({0:{0:QQ(1), 2:QQ(-1)}, 1:{1:QQ(1), 2:QQ(2)}}, (3, 3), QQ)
    assert A.rref() == (Arref, [0, 1])

    A = SDM({0:{0:QQ(1), 1:QQ(2), 3:QQ(1)},
             1:{0:QQ(1), 1:QQ(1), 2:QQ(9)}}, (2, 4), QQ)
    Arref = SDM({0:{0:QQ(1), 2:QQ(18), 3:QQ(-1)},
                 1:{1:QQ(1), 2:QQ(-9), 3:QQ(1)}}, (2, 4), QQ)
    assert A.rref() == (Arref, [0, 1])

    A = SDM({0:{0:QQ(1), 1:QQ(1), 2:QQ(1)},
             1:{0:QQ(1), 1:QQ(2), 2:QQ(2)}}, (2, 3), QQ)
    Arref = SDM(
            {0: {0: QQ(1,1)}, 1: {1: QQ(1,1), 2: QQ(1,1)}},
            (2, 3), QQ)
    assert A.rref() == (Arref, [0, 1])
