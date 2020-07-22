from sympy.core.compatibility import HAS_GMPY
from sympy.core.numbers import Rational
from sympy.matrices.common import (NonInvertibleMatrixError,
    NonSquareMatrixError, ShapeError)
from sympy.matrices.dense import Matrix
from sympy.polys import ZZ, QQ
from sympy.polys.domainmatrix import (
        DDM,
        DDMBadInputError, DDMDomainError, DDMShapeError,
        DomainMatrix,
        ddm_iadd, ddm_isub, ddm_ineg, ddm_imatmul,
        ddm_irref, ddm_idet, ddm_iinv, ddm_ilu, ddm_ilu_split, ddm_ilu_solve,
        ddm_berk,
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


def test_ddm_iadd():
    a = [[1, 2], [3, 4]]
    b = [[5, 6], [7, 8]]
    ddm_iadd(a, b)
    assert a == [[6, 8], [10, 12]]


def test_ddm_isub():
    a = [[1, 2], [3, 4]]
    b = [[5, 6], [7, 8]]
    ddm_isub(a, b)
    assert a == [[-4, -4], [-4, -4]]


def test_ddm_ineg():
    a = [[1, 2], [3, 4]]
    ddm_ineg(a)
    assert a == [[-1, -2], [-3, -4]]


def test_ddm_imatmul():
    a = [[1, 2, 3], [4, 5, 6]]
    b = [[1, 2], [3, 4], [5, 6]]

    c1 = [[0, 0], [0, 0]]
    ddm_imatmul(c1, a, b)
    assert c1 == [[22, 28], [49, 64]]

    c2 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    ddm_imatmul(c2, b, a)
    assert c2 == [[9, 12, 15], [19, 26, 33], [29, 40, 51]]

    b3 = [[1], [2], [3]]
    c3 = [[0], [0]]
    ddm_imatmul(c3, a, b3)
    assert c3 == [[14], [32]]


def test_ddm_irref():
    # Empty matrix
    A = []
    Ar = []
    pivots = []
    assert ddm_irref(A) == pivots
    assert A == Ar

    # Standard square case
    A = [[QQ(0), QQ(1)], [QQ(1), QQ(1)]]
    Ar = [[QQ(1), QQ(0)], [QQ(0), QQ(1)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    # m < n  case
    A = [[QQ(1), QQ(2), QQ(1)], [QQ(3), QQ(4), QQ(1)]]
    Ar = [[QQ(1), QQ(0), QQ(-1)], [QQ(0), QQ(1), QQ(1)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    # same m < n  but reversed
    A = [[QQ(3), QQ(4), QQ(1)], [QQ(1), QQ(2), QQ(1)]]
    Ar = [[QQ(1), QQ(0), QQ(-1)], [QQ(0), QQ(1), QQ(1)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    # m > n case
    A = [[QQ(1), QQ(0)], [QQ(1), QQ(3)], [QQ(0), QQ(1)]]
    Ar = [[QQ(1), QQ(0)], [QQ(0), QQ(1)], [QQ(0), QQ(0)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    # Example with missing pivot
    A = [[QQ(1), QQ(0), QQ(1)], [QQ(3), QQ(0), QQ(1)]]
    Ar = [[QQ(1), QQ(0), QQ(0)], [QQ(0), QQ(0), QQ(1)]]
    pivots = [0, 2]
    assert ddm_irref(A) == pivots
    assert A == Ar

    # Example with missing pivot and no replacement
    A = [[QQ(0), QQ(1)], [QQ(0), QQ(2)], [QQ(1), QQ(0)]]
    Ar = [[QQ(1), QQ(0)], [QQ(0), QQ(1)], [QQ(0), QQ(0)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar


def test_ddm_idet():
    A = []
    assert ddm_idet(A, ZZ) == ZZ(1)

    A = [[ZZ(2)]]
    assert ddm_idet(A, ZZ) == ZZ(2)

    A = [[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]]
    assert ddm_idet(A, ZZ) == ZZ(-2)

    A = [[ZZ(1), ZZ(2), ZZ(3)], [ZZ(1), ZZ(2), ZZ(4)], [ZZ(1), ZZ(3), ZZ(5)]]
    assert ddm_idet(A, ZZ) == ZZ(-1)

    A = [[ZZ(1), ZZ(2), ZZ(3)], [ZZ(1), ZZ(2), ZZ(4)], [ZZ(1), ZZ(2), ZZ(5)]]
    assert ddm_idet(A, ZZ) == ZZ(0)

    A = [[QQ(1, 2), QQ(1, 2)], [QQ(1, 3), QQ(1, 4)]]
    assert ddm_idet(A, QQ) == QQ(-1, 24)


def test_ddm_inv():
    A = []
    Ainv = []
    ddm_iinv(Ainv, A, QQ)
    assert Ainv == A

    A = []
    Ainv = []
    raises(ValueError, lambda: ddm_iinv(Ainv, A, ZZ))

    A = [[QQ(1), QQ(2)]]
    Ainv = [[QQ(0), QQ(0)]]
    raises(NonSquareMatrixError, lambda: ddm_iinv(Ainv, A, QQ))

    A = [[QQ(1, 1), QQ(2, 1)], [QQ(3, 1), QQ(4, 1)]]
    Ainv = [[QQ(0), QQ(0)], [QQ(0), QQ(0)]]
    Ainv_expected = [[QQ(-2, 1), QQ(1, 1)], [QQ(3, 2), QQ(-1, 2)]]
    ddm_iinv(Ainv, A, QQ)
    assert Ainv == Ainv_expected

    A = [[QQ(1, 1), QQ(2, 1)], [QQ(2, 1), QQ(4, 1)]]
    Ainv = [[QQ(0), QQ(0)], [QQ(0), QQ(0)]]
    raises(NonInvertibleMatrixError, lambda: ddm_iinv(Ainv, A, QQ))


def test_ddm_ilu():
    A = []
    Alu = []
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == []

    A = [[]]
    Alu = [[]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == []

    A = [[QQ(1), QQ(2)], [QQ(3), QQ(4)]]
    Alu = [[QQ(1), QQ(2)], [QQ(3), QQ(-2)]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == []

    A = [[QQ(0), QQ(2)], [QQ(3), QQ(4)]]
    Alu = [[QQ(3), QQ(4)], [QQ(0), QQ(2)]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == [(0, 1)]

    A = [[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(5), QQ(6)], [QQ(7), QQ(8), QQ(9)]]
    Alu = [[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(-3), QQ(-6)], [QQ(7), QQ(2), QQ(0)]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == []

    A = [[QQ(0), QQ(1), QQ(2)], [QQ(0), QQ(1), QQ(3)], [QQ(1), QQ(1), QQ(2)]]
    Alu = [[QQ(1), QQ(1), QQ(2)], [QQ(0), QQ(1), QQ(3)], [QQ(0), QQ(1), QQ(-1)]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == [(0, 2)]

    A = [[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(5), QQ(6)]]
    Alu = [[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(-3), QQ(-6)]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == []

    A = [[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]]
    Alu = [[QQ(1), QQ(2)], [QQ(3), QQ(-2)], [QQ(5), QQ(2)]]
    swaps = ddm_ilu(A)
    assert A == Alu
    assert swaps == []


def test_ddm_ilu_split():
    U = []
    L = []
    Uexp = []
    Lexp = []
    swaps = ddm_ilu_split(L, U, QQ)
    assert U == Uexp
    assert L == Lexp
    assert swaps == []

    U = [[]]
    L = [[QQ(1)]]
    Uexp = [[]]
    Lexp = [[QQ(1)]]
    swaps = ddm_ilu_split(L, U, QQ)
    assert U == Uexp
    assert L == Lexp
    assert swaps == []

    U = [[QQ(1), QQ(2)], [QQ(3), QQ(4)]]
    L = [[QQ(1), QQ(0)], [QQ(0), QQ(1)]]
    Uexp = [[QQ(1), QQ(2)], [QQ(0), QQ(-2)]]
    Lexp = [[QQ(1), QQ(0)], [QQ(3), QQ(1)]]
    swaps = ddm_ilu_split(L, U, QQ)
    assert U == Uexp
    assert L == Lexp
    assert swaps == []

    U = [[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(5), QQ(6)]]
    L = [[QQ(1), QQ(0)], [QQ(0), QQ(1)]]
    Uexp = [[QQ(1), QQ(2), QQ(3)], [QQ(0), QQ(-3), QQ(-6)]]
    Lexp = [[QQ(1), QQ(0)], [QQ(4), QQ(1)]]
    swaps = ddm_ilu_split(L, U, QQ)
    assert U == Uexp
    assert L == Lexp
    assert swaps == []

    U = [[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]]
    L = [[QQ(1), QQ(0), QQ(0)], [QQ(0), QQ(1), QQ(0)], [QQ(0), QQ(0), QQ(1)]]
    Uexp = [[QQ(1), QQ(2)], [QQ(0), QQ(-2)], [QQ(0), QQ(0)]]
    Lexp = [[QQ(1), QQ(0), QQ(0)], [QQ(3), QQ(1), QQ(0)], [QQ(5), QQ(2), QQ(1)]]
    swaps = ddm_ilu_split(L, U, QQ)
    assert U == Uexp
    assert L == Lexp
    assert swaps == []


def test_ddm_ilu_solve():
    # Basic example
    # A = [[QQ(1), QQ(2)], [QQ(3), QQ(4)]]
    U = [[QQ(1), QQ(2)], [QQ(0), QQ(-2)]]
    L = [[QQ(1), QQ(0)], [QQ(3), QQ(1)]]
    swaps = []
    b = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    x = DDM([[QQ(0)], [QQ(0)]], (2, 1), QQ)
    xexp = DDM([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    ddm_ilu_solve(x, L, U, swaps, b)
    assert x == xexp

    # Example with swaps
    # A = [[QQ(0), QQ(2)], [QQ(3), QQ(4)]]
    U = [[QQ(3), QQ(4)], [QQ(0), QQ(2)]]
    L = [[QQ(1), QQ(0)], [QQ(0), QQ(1)]]
    swaps = [(0, 1)]
    b = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    x = DDM([[QQ(0)], [QQ(0)]], (2, 1), QQ)
    xexp = DDM([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    ddm_ilu_solve(x, L, U, swaps, b)
    assert x == xexp

    # Overdetermined, consistent
    # A = DDM([[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]], (3, 2), QQ)
    U = [[QQ(1), QQ(2)], [QQ(0), QQ(-2)], [QQ(0), QQ(0)]]
    L = [[QQ(1), QQ(0), QQ(0)], [QQ(3), QQ(1), QQ(0)], [QQ(5), QQ(2), QQ(1)]]
    swaps = []
    b = DDM([[QQ(1)], [QQ(2)], [QQ(3)]], (3, 1), QQ)
    x = DDM([[QQ(0)], [QQ(0)]], (2, 1), QQ)
    xexp = DDM([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    ddm_ilu_solve(x, L, U, swaps, b)
    assert x == xexp

    # Overdetermined, inconsistent
    b = DDM([[QQ(1)], [QQ(2)], [QQ(4)]], (3, 1), QQ)
    raises(NonInvertibleMatrixError, lambda: ddm_ilu_solve(x, L, U, swaps, b))

    # Square, noninvertible
    # A = DDM([[QQ(1), QQ(2)], [QQ(1), QQ(2)]], (2, 2), QQ)
    U = [[QQ(1), QQ(2)], [QQ(0), QQ(0)]]
    L = [[QQ(1), QQ(0)], [QQ(1), QQ(1)]]
    swaps = []
    b = DDM([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    raises(NonInvertibleMatrixError, lambda: ddm_ilu_solve(x, L, U, swaps, b))

    # Underdetermined
    # A = DDM([[QQ(1), QQ(2)]], (1, 2), QQ)
    U = [[QQ(1), QQ(2)]]
    L = [[QQ(1)]]
    swaps = []
    b = DDM([[QQ(3)]], (1, 1), QQ)
    raises(NotImplementedError, lambda: ddm_ilu_solve(x, L, U, swaps, b))

    # Shape mismatch
    b3 = DDM([[QQ(1)], [QQ(2)], [QQ(3)]], (3, 1), QQ)
    raises(DDMShapeError, lambda: ddm_ilu_solve(x, L, U, swaps, b3))

    # Empty shape mismatch
    U = [[QQ(1)]]
    L = [[QQ(1)]]
    swaps = []
    x = [[QQ(1)]]
    b = []
    raises(DDMShapeError, lambda: ddm_ilu_solve(x, L, U, swaps, b))

    # Empty system
    U = []
    L = []
    swaps = []
    b = []
    x = []
    ddm_ilu_solve(x, L, U, swaps, b)
    assert x == []


def test_ddm_charpoly():
    A = []
    assert ddm_berk(A, ZZ) == [[ZZ(1)]]

    A = [[ZZ(1), ZZ(2), ZZ(3)], [ZZ(4), ZZ(5), ZZ(6)], [ZZ(7), ZZ(8), ZZ(9)]]
    Avec = [[ZZ(1)], [ZZ(-15)], [ZZ(-18)], [ZZ(0)]]
    assert ddm_berk(A, ZZ) == Avec

    A = DDM([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
    raises(DDMShapeError, lambda: ddm_berk(A, ZZ))


def test_DomainMatrix_init():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.rep == DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    raises(DDMBadInputError, lambda: DomainMatrix([[ZZ(1), ZZ(2)]], (2, 2), ZZ))


def test_DomainMatrix_from_ddm():
    ddm = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A = DomainMatrix.from_ddm(ddm)
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == ZZ


def test_DomainMatrix_from_list_sympy():
    ddm = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A = DomainMatrix.from_list_sympy(2, 2, [[1, 2], [3, 4]])
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == ZZ


def test_DomainMatrix_eq():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A == A
    B = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(1)]], (2, 2), ZZ)
    assert A != B
    C = [[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]]
    assert A != C


def test_DomainMatrix_get_domain():
    K, items = DomainMatrix.get_domain([1, 2, 3, 4])
    assert items == [ZZ(1), ZZ(2), ZZ(3), ZZ(4)]
    assert K == ZZ

    K, items = DomainMatrix.get_domain([1, 2, 3, Rational(1, 2)])
    assert items == [QQ(1), QQ(2), QQ(3), QQ(1, 2)]
    assert K == QQ


def test_DomainMatrix_convert_to():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = A.convert_to(QQ)
    assert Aq == DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)


def test_DomainMatrix_to_field():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = A.to_field()
    assert Aq == DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)


def test_DomainMatrix_unify():
    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    assert Az.unify(Az) == (Az, Az)
    assert Az.unify(Aq) == (Aq, Aq)
    assert Aq.unify(Az) == (Aq, Aq)
    assert Aq.unify(Aq) == (Aq, Aq)


def test_DomainMatrix_to_Matrix():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.to_Matrix() == Matrix([[1, 2], [3, 4]])


def test_DomainMatrix_repr():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert repr(A) == 'DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)'


def test_DomainMatrix_add():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    B = DomainMatrix([[ZZ(2), ZZ(4)], [ZZ(6), ZZ(8)]], (2, 2), ZZ)
    assert A + A == A.add(A) == B

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    L = [[2, 3], [3, 4]]
    raises(TypeError, lambda: A + L)
    raises(TypeError, lambda: L + A)

    A1 = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A2 = DomainMatrix([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
    raises(ShapeError, lambda: A1 + A2)
    raises(ShapeError, lambda: A2 + A1)
    raises(ShapeError, lambda: A1.add(A2))
    raises(ShapeError, lambda: A2.add(A1))

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    raises(ValueError, lambda: Az + Aq)
    raises(ValueError, lambda: Aq + Az)
    raises(ValueError, lambda: Az.add(Aq))
    raises(ValueError, lambda: Aq.add(Az))


def test_DomainMatrix_sub():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    B = DomainMatrix([[ZZ(0), ZZ(0)], [ZZ(0), ZZ(0)]], (2, 2), ZZ)
    assert A - A == A.sub(A) == B

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    L = [[2, 3], [3, 4]]
    raises(TypeError, lambda: A - L)
    raises(TypeError, lambda: L - A)

    A1 = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A2 = DomainMatrix([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
    raises(ShapeError, lambda: A1 - A2)
    raises(ShapeError, lambda: A2 - A1)
    raises(ShapeError, lambda: A1.sub(A2))
    raises(ShapeError, lambda: A2.sub(A1))

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    raises(ValueError, lambda: Az - Aq)
    raises(ValueError, lambda: Aq - Az)
    raises(ValueError, lambda: Az.sub(Aq))
    raises(ValueError, lambda: Aq.sub(Az))


def test_DomainMatrix_neg():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aneg = DomainMatrix([[ZZ(-1), ZZ(-2)], [ZZ(-3), ZZ(-4)]], (2, 2), ZZ)
    assert -A == A.neg() == Aneg


def test_DomainMatrix_mul():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A2 = DomainMatrix([[ZZ(7), ZZ(10)], [ZZ(15), ZZ(22)]], (2, 2), ZZ)
    assert A*A == A.matmul(A) == A2

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    L = [[1, 2], [3, 4]]
    raises(TypeError, lambda: A * L)
    raises(TypeError, lambda: L * A)

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    raises(DDMDomainError, lambda: Az * Aq)
    raises(DDMDomainError, lambda: Aq * Az)
    raises(DDMDomainError, lambda: Az.matmul(Aq))
    raises(DDMDomainError, lambda: Aq.matmul(Az))


def test_DomainMatrix_pow():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A2 = DomainMatrix([[ZZ(7), ZZ(10)], [ZZ(15), ZZ(22)]], (2, 2), ZZ)
    A3 = DomainMatrix([[ZZ(37), ZZ(54)], [ZZ(81), ZZ(118)]], (2, 2), ZZ)
    eye = DomainMatrix([[ZZ(1), ZZ(0)], [ZZ(0), ZZ(1)]], (2, 2), ZZ)
    assert A**0 == A.pow(0) == eye
    assert A**1 == A.pow(1) == A
    assert A**2 == A.pow(2) == A2
    assert A**3 == A.pow(3) == A3

    raises(TypeError, lambda: A ** Rational(1, 2))
    raises(NotImplementedError, lambda: A ** -1)
    raises(NotImplementedError, lambda: A.pow(-1))


def test_DomainMatrix_rref():
    A = DomainMatrix([], (0, 1), QQ)
    assert A.rref() == (A, ())

    A = DomainMatrix([[QQ(1)]], (1, 1), QQ)
    assert A.rref() == (A, (0,))

    A = DomainMatrix([[QQ(0)]], (1, 1), QQ)
    assert A.rref() == (A, ())

    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    Ar, pivots = A.rref()
    assert Ar == DomainMatrix([[QQ(1), QQ(0)], [QQ(0), QQ(1)]], (2, 2), QQ)
    assert pivots == (0, 1)

    A = DomainMatrix([[QQ(0), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    Ar, pivots = A.rref()
    assert Ar == DomainMatrix([[QQ(1), QQ(0)], [QQ(0), QQ(1)]], (2, 2), QQ)
    assert pivots == (0, 1)

    A = DomainMatrix([[QQ(0), QQ(2)], [QQ(0), QQ(4)]], (2, 2), QQ)
    Ar, pivots = A.rref()
    assert Ar == DomainMatrix([[QQ(0), QQ(1)], [QQ(0), QQ(0)]], (2, 2), QQ)
    assert pivots == (1,)

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    raises(ValueError, lambda: Az.rref())


def test_DomainMatrix_inv():
    A = DomainMatrix([], (0, 0), QQ)
    assert A.inv() == A

    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    Ainv = DomainMatrix([[QQ(-2), QQ(1)], [QQ(3, 2), QQ(-1, 2)]], (2, 2), QQ)
    assert A.inv() == Ainv

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    raises(ValueError, lambda: Az.inv())

    Ans = DomainMatrix([[QQ(1), QQ(2)]], (1, 2), QQ)
    raises(NonSquareMatrixError, lambda: Ans.inv())

    Aninv = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(6)]], (2, 2), QQ)
    raises(NonInvertibleMatrixError, lambda: Aninv.inv())


def test_DomainMatrix_det():
    A = DomainMatrix([], (0, 0), ZZ)
    assert A.det() == 1

    A = DomainMatrix([[1]], (1, 1), ZZ)
    assert A.det() == 1

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.det() == ZZ(-2)

    A = DomainMatrix([[ZZ(1), ZZ(2), ZZ(3)], [ZZ(1), ZZ(2), ZZ(4)], [ZZ(1), ZZ(3), ZZ(5)]], (3, 3), ZZ)
    assert A.det() == ZZ(-1)

    A = DomainMatrix([[ZZ(1), ZZ(2), ZZ(3)], [ZZ(1), ZZ(2), ZZ(4)], [ZZ(1), ZZ(2), ZZ(5)]], (3, 3), ZZ)
    assert A.det() == ZZ(0)

    Ans = DomainMatrix([[QQ(1), QQ(2)]], (1, 2), QQ)
    raises(NonSquareMatrixError, lambda: Ans.det())

    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    assert A.det() == QQ(-2)


def test_DomainMatrix_lu():
    A = DomainMatrix([], (0, 0), QQ)
    assert A.lu() == (A, A, [])

    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    L = DomainMatrix([[QQ(1), QQ(0)], [QQ(3), QQ(1)]], (2, 2), QQ)
    U = DomainMatrix([[QQ(1), QQ(2)], [QQ(0), QQ(-2)]], (2, 2), QQ)
    swaps = []
    assert A.lu() == (L, U, swaps)

    A = DomainMatrix([[QQ(0), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    L = DomainMatrix([[QQ(1), QQ(0)], [QQ(0), QQ(1)]], (2, 2), QQ)
    U = DomainMatrix([[QQ(3), QQ(4)], [QQ(0), QQ(2)]], (2, 2), QQ)
    swaps = [(0, 1)]
    assert A.lu() == (L, U, swaps)

    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(2), QQ(4)]], (2, 2), QQ)
    L = DomainMatrix([[QQ(1), QQ(0)], [QQ(2), QQ(1)]], (2, 2), QQ)
    U = DomainMatrix([[QQ(1), QQ(2)], [QQ(0), QQ(0)]], (2, 2), QQ)
    swaps = []
    assert A.lu() == (L, U, swaps)

    A = DomainMatrix([[QQ(0), QQ(2)], [QQ(0), QQ(4)]], (2, 2), QQ)
    L = DomainMatrix([[QQ(1), QQ(0)], [QQ(0), QQ(1)]], (2, 2), QQ)
    U = DomainMatrix([[QQ(0), QQ(2)], [QQ(0), QQ(4)]], (2, 2), QQ)
    swaps = []
    assert A.lu() == (L, U, swaps)

    A = DomainMatrix([[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(5), QQ(6)]], (2, 3), QQ)
    L = DomainMatrix([[QQ(1), QQ(0)], [QQ(4), QQ(1)]], (2, 2), QQ)
    U = DomainMatrix([[QQ(1), QQ(2), QQ(3)], [QQ(0), QQ(-3), QQ(-6)]], (2, 3), QQ)
    swaps = []
    assert A.lu() == (L, U, swaps)

    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]], (3, 2), QQ)
    L = DomainMatrix([
        [QQ(1), QQ(0), QQ(0)],
        [QQ(3), QQ(1), QQ(0)],
        [QQ(5), QQ(2), QQ(1)]], (3, 3), QQ)
    U = DomainMatrix([[QQ(1), QQ(2)], [QQ(0), QQ(-2)], [QQ(0), QQ(0)]], (3, 2), QQ)
    swaps = []
    assert A.lu() == (L, U, swaps)

    A = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]]
    L = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 1]]
    U = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]]
    to_dom = lambda rows, dom: [[dom(e) for e in row] for row in rows]
    A = DomainMatrix(to_dom(A, QQ), (4, 4), QQ)
    L = DomainMatrix(to_dom(L, QQ), (4, 4), QQ)
    U = DomainMatrix(to_dom(U, QQ), (4, 4), QQ)
    assert A.lu() == (L, U, [])

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    raises(ValueError, lambda: A.lu())


def test_DomainMatrix_lu_solve():
    # Base case
    A = b = x = DomainMatrix([], (0, 0), QQ)
    assert A.lu_solve(b) == x

    # Basic example
    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    b = DomainMatrix([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    x = DomainMatrix([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    assert A.lu_solve(b) == x

    # Example with swaps
    A = DomainMatrix([[QQ(0), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    b = DomainMatrix([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    x = DomainMatrix([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    assert A.lu_solve(b) == x

    # Non-invertible
    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(2), QQ(4)]], (2, 2), QQ)
    b = DomainMatrix([[QQ(1)], [QQ(2)]], (2, 1), QQ)
    raises(NonInvertibleMatrixError, lambda: A.lu_solve(b))

    # Overdetermined, consistent
    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]], (3, 2), QQ)
    b = DomainMatrix([[QQ(1)], [QQ(2)], [QQ(3)]], (3, 1), QQ)
    x = DomainMatrix([[QQ(0)], [QQ(1, 2)]], (2, 1), QQ)
    assert A.lu_solve(b) == x

    # Overdetermined, inconsistent
    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)], [QQ(5), QQ(6)]], (3, 2), QQ)
    b = DomainMatrix([[QQ(1)], [QQ(2)], [QQ(4)]], (3, 1), QQ)
    raises(NonInvertibleMatrixError, lambda: A.lu_solve(b))

    # Underdetermined
    A = DomainMatrix([[QQ(1), QQ(2)]], (1, 2), QQ)
    b = DomainMatrix([[QQ(1)]], (1, 1), QQ)
    raises(NotImplementedError, lambda: A.lu_solve(b))

    # Non-field
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    b = DomainMatrix([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    raises(ValueError, lambda: A.lu_solve(b))

    # Shape mismatch
    A = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    b = DomainMatrix([[QQ(1), QQ(2)]], (1, 2), QQ)
    raises(ShapeError, lambda: A.lu_solve(b))


def test_DomainMatrix_charpoly():
    A = DomainMatrix([], (0, 0), ZZ)
    assert A.charpoly() == [ZZ(1)]

    A = DomainMatrix([[1]], (1, 1), ZZ)
    assert A.charpoly() == [ZZ(1), ZZ(-1)]

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.charpoly() == [ZZ(1), ZZ(-5), ZZ(-2)]

    A = DomainMatrix([[ZZ(1), ZZ(2), ZZ(3)], [ZZ(4), ZZ(5), ZZ(6)], [ZZ(7), ZZ(8), ZZ(9)]], (3, 3), ZZ)
    assert A.charpoly() == [ZZ(1), ZZ(-15), ZZ(-18), ZZ(0)]

    Ans = DomainMatrix([[QQ(1), QQ(2)]], (1, 2), QQ)
    raises(NonSquareMatrixError, lambda: Ans.charpoly())
