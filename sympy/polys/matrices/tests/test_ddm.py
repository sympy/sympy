from sympy.testing.pytest import raises
from sympy.core.compatibility import HAS_GMPY

from sympy.polys import ZZ, QQ

from sympy.polys.matrices.ddm import DDM
from sympy.polys.matrices.exceptions import (
        DDMShapeError, NonInvertibleMatrixError, DDMDomainError,
        DDMBadInputError)


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


def test_DDM_getsetitem():
    ddm = DDM([[ZZ(2), ZZ(3)], [ZZ(4), ZZ(5)]], (2, 2), ZZ)

    assert ddm[0][0] == ZZ(2)
    assert ddm[0][1] == ZZ(3)
    assert ddm[1][0] == ZZ(4)
    assert ddm[1][1] == ZZ(5)

    raises(IndexError, lambda: ddm[2][0])
    raises(IndexError, lambda: ddm[0][2])

    ddm[0][0] = ZZ(-1)
    assert ddm[0][0] == ZZ(-1)


def test_DDM_str():
    ddm = DDM([[ZZ(0), ZZ(1)], [ZZ(2), ZZ(3)]], (2, 2), ZZ)
    if HAS_GMPY: # pragma: no cover
        assert str(ddm) == '[[0, 1], [2, 3]]'
        assert repr(ddm) == 'DDM([[mpz(0), mpz(1)], [mpz(2), mpz(3)]], (2, 2), ZZ)'
    else:        # pragma: no cover
        assert repr(ddm) == 'DDM([[0, 1], [2, 3]], (2, 2), ZZ)'
        assert str(ddm) == '[[0, 1], [2, 3]]'


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


def test_DDM_convert_to():
    ddm = DDM([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
    assert ddm.convert_to(ZZ) == ddm
    ddmq = ddm.convert_to(QQ)
    assert ddmq.domain == QQ


def test_DDM_zeros():
    ddmz = DDM.zeros((3, 4), QQ)
    assert list(ddmz) == [[QQ(0)] * 4] * 3
    assert ddmz.shape == (3, 4)
    assert ddmz.domain == QQ


def test_DDM_eye():
    ddmz = DDM.eye(3, QQ)
    f = lambda i, j: QQ(1) if i == j else QQ(0)
    assert list(ddmz) == [[f(i, j) for i in range(3)] for j in range(3)]
    assert ddmz.shape == (3, 3)
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
    A2 = DDM([[ZZ(2)]], (1, 1), ZZ)
    assert A * ZZ(2) == A2
    assert ZZ(2) * A == A2
    raises(TypeError, lambda: [[1]] * A)
    raises(TypeError, lambda: A * [[1]])


def test_DDM_matmul():
    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    B = DDM([[ZZ(3), ZZ(4)]], (1, 2), ZZ)
    AB = DDM([[ZZ(3), ZZ(4)], [ZZ(6), ZZ(8)]], (2, 2), ZZ)
    BA = DDM([[ZZ(11)]], (1, 1), ZZ)

    assert A @ B == A.matmul(B) == AB
    assert B @ A == B.matmul(A) == BA

    raises(TypeError, lambda: A @ 1)
    raises(TypeError, lambda: A @ [[3, 4]])

    Bq = DDM([[QQ(3), QQ(4)]], (1, 2), QQ)

    raises(DDMDomainError, lambda: A @ Bq)
    raises(DDMDomainError, lambda: Bq @ A)

    C = DDM([[ZZ(1)]], (1, 1), ZZ)

    assert A @ C == A.matmul(C) == A

    raises(DDMShapeError, lambda: C @ A)
    raises(DDMShapeError, lambda: C.matmul(A))

    Z04 = DDM([], (0, 4), ZZ)
    Z40 = DDM([[]]*4, (4, 0), ZZ)
    Z50 = DDM([[]]*5, (5, 0), ZZ)
    Z05 = DDM([], (0, 5), ZZ)
    Z45 = DDM([[0] * 5] * 4, (4, 5), ZZ)
    Z54 = DDM([[0] * 4] * 5, (5, 4), ZZ)
    Z00 = DDM([], (0, 0), ZZ)

    assert Z04 @ Z45 == Z04.matmul(Z45) == Z05
    assert Z45 @ Z50 == Z45.matmul(Z50) == Z40
    assert Z00 @ Z04 == Z00.matmul(Z04) == Z04
    assert Z50 @ Z00 == Z50.matmul(Z00) == Z50
    assert Z00 @ Z00 == Z00.matmul(Z00) == Z00
    assert Z50 @ Z04 == Z50.matmul(Z04) == Z54

    raises(DDMShapeError, lambda: Z05 @ Z40)
    raises(DDMShapeError, lambda: Z05.matmul(Z40))

def test_DDM_hstack():

    A = DDM([[ZZ(1), ZZ(2), ZZ(3)]], (1, 3), ZZ)
    B = DDM([[ZZ(4), ZZ(5)]], (1, 2), ZZ)
    Ah = A.hstack(B)

    assert Ah.shape == (1, 5)
    assert Ah.domain == ZZ
    assert Ah == DDM([[ZZ(1), ZZ(2), ZZ(3), ZZ(4), ZZ(5)]], (1, 5), ZZ)

def test_DDM_rref():

    A = DDM([], (0, 4), QQ)
    assert A.rref() == (A, [])

    A = DDM([[QQ(0), QQ(1)], [QQ(1), QQ(1)]], (2, 2), QQ)
    Ar = DDM([[QQ(1), QQ(0)], [QQ(0), QQ(1)]], (2, 2), QQ)
    pivots = [0, 1]
    assert A.rref() == (Ar, pivots)

    A = DDM([[QQ(1), QQ(2), QQ(1)], [QQ(3), QQ(4), QQ(1)]], (2, 3), QQ)
    Ar = DDM([[QQ(1), QQ(0), QQ(-1)], [QQ(0), QQ(1), QQ(1)]], (2, 3), QQ)
    pivots = [0, 1]
    assert A.rref() == (Ar, pivots)

    A = DDM([[QQ(3), QQ(4), QQ(1)], [QQ(1), QQ(2), QQ(1)]], (2, 3), QQ)
    Ar = DDM([[QQ(1), QQ(0), QQ(-1)], [QQ(0), QQ(1), QQ(1)]], (2, 3), QQ)
    pivots = [0, 1]
    assert A.rref() == (Ar, pivots)

    A = DDM([[QQ(1), QQ(0)], [QQ(1), QQ(3)], [QQ(0), QQ(1)]], (3, 2), QQ)
    Ar = DDM([[QQ(1), QQ(0)], [QQ(0), QQ(1)], [QQ(0), QQ(0)]], (3, 2), QQ)
    pivots = [0, 1]
    assert A.rref() == (Ar, pivots)

    A = DDM([[QQ(1), QQ(0), QQ(1)], [QQ(3), QQ(0), QQ(1)]], (2, 3), QQ)
    Ar = DDM([[QQ(1), QQ(0), QQ(0)], [QQ(0), QQ(0), QQ(1)]], (2, 3), QQ)
    pivots = [0, 2]
    assert A.rref() == (Ar, pivots)


def test_DDM_nullspace():
     A = DDM([[QQ(1), QQ(1)], [QQ(1), QQ(1)]], (2, 2), QQ)
     Anull = DDM([[QQ(-1), QQ(1)]], (1, 2), QQ)
     nonpivots = [1]
     assert A.nullspace() == (Anull, nonpivots)


def test_DDM_det():
    # 0x0 case
    A = DDM([], (0, 0), ZZ)
    assert A.det() == ZZ(1)

    # 1x1 case
    A = DDM([[ZZ(2)]], (1, 1), ZZ)
    assert A.det() == ZZ(2)

    # 2x2 case
    A = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.det() == ZZ(-2)

    # 3x3 with swap
    A = DDM([[ZZ(1), ZZ(2), ZZ(3)], [ZZ(1), ZZ(2), ZZ(4)], [ZZ(1), ZZ(2), ZZ(5)]], (3, 3), ZZ)
    assert A.det() == ZZ(0)

    # 2x2 QQ case
    A = DDM([[QQ(1, 2), QQ(1, 2)], [QQ(1, 3), QQ(1, 4)]], (2, 2), QQ)
    assert A.det() == QQ(-1, 24)

    # Nonsquare error
    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    raises(DDMShapeError, lambda: A.det())

    # Nonsquare error with empty matrix
    A = DDM([], (0, 1), ZZ)
    raises(DDMShapeError, lambda: A.det())


def test_DDM_inv():
    A = DDM([[QQ(1, 1), QQ(2, 1)], [QQ(3, 1), QQ(4, 1)]], (2, 2), QQ)
    Ainv = DDM([[QQ(-2, 1), QQ(1, 1)], [QQ(3, 2), QQ(-1, 2)]], (2, 2), QQ)
    assert A.inv() == Ainv

    A = DDM([[QQ(1), QQ(2)]], (1, 2), QQ)
    raises(DDMShapeError, lambda: A.inv())

    A = DDM([[ZZ(2)]], (1, 1), ZZ)
    raises(ValueError, lambda: A.inv())

    A = DDM([], (0, 0), QQ)
    assert A.inv() == A

    A = DDM([[QQ(1), QQ(2)], [QQ(2), QQ(4)]], (2, 2), QQ)
    raises(NonInvertibleMatrixError, lambda: A.inv())


def test_DDM_lu():
    A = DDM([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    L, U, swaps = A.lu()
    assert L == DDM([[QQ(1), QQ(0)], [QQ(3), QQ(1)]], (2, 2), QQ)
    assert U == DDM([[QQ(1), QQ(2)], [QQ(0), QQ(-2)]], (2, 2), QQ)
    assert swaps == []

    A = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]]
    Lexp = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 1]]
    Uexp = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]]
    to_dom = lambda rows, dom: [[dom(e) for e in row] for row in rows]
    A = DDM(to_dom(A, QQ), (4, 4), QQ)
    Lexp = DDM(to_dom(Lexp, QQ), (4, 4), QQ)
    Uexp = DDM(to_dom(Uexp, QQ), (4, 4), QQ)
    L, U, swaps = A.lu()
    assert L == Lexp
    assert U == Uexp
    assert swaps == []


def test_DDM_lu_solve():
    # Basic example
    A = DDM([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    b = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    x = DDM([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    assert A.lu_solve(b) == x

    # Example with swaps
    A = DDM([[QQ(0), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    assert A.lu_solve(b) == x

    # Overdetermined, consistent
    A = DDM([[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]], (3, 2), QQ)
    b = DDM([[QQ(1)], [QQ(2)], [QQ(3)]], (3, 1), QQ)
    assert A.lu_solve(b) == x

    # Overdetermined, inconsistent
    b = DDM([[QQ(1)], [QQ(2)], [QQ(4)]], (3, 1), QQ)
    raises(NonInvertibleMatrixError, lambda: A.lu_solve(b))

    # Square, noninvertible
    A = DDM([[QQ(1), QQ(2)], [QQ(1), QQ(2)]], (2, 2), QQ)
    b = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    raises(NonInvertibleMatrixError, lambda: A.lu_solve(b))

    # Underdetermined
    A = DDM([[QQ(1), QQ(2)]], (1, 2), QQ)
    b = DDM([[QQ(3)]], (1, 1), QQ)
    raises(NotImplementedError, lambda: A.lu_solve(b))

    # Domain mismatch
    bz = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    raises(DDMDomainError, lambda: A.lu_solve(bz))

    # Shape mismatch
    b3 = DDM([[QQ(1)], [QQ(2)], [QQ(3)]], (3, 1), QQ)
    raises(DDMShapeError, lambda: A.lu_solve(b3))


def test_DDM_charpoly():
    A = DDM([], (0, 0), ZZ)
    assert A.charpoly() == [ZZ(1)]

    A = DDM([
        [ZZ(1), ZZ(2), ZZ(3)],
        [ZZ(4), ZZ(5), ZZ(6)],
        [ZZ(7), ZZ(8), ZZ(9)]], (3, 3), ZZ)
    Avec = [ZZ(1), ZZ(-15), ZZ(-18), ZZ(0)]
    assert A.charpoly() == Avec

    A = DDM([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
    raises(DDMShapeError, lambda: A.charpoly())
