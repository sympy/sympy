from sympy.core.compatibility import HAS_GMPY
from sympy.matrices.common import NonInvertibleMatrixError, NonSquareMatrixError
from sympy.polys import ZZ, QQ
from sympy.polys.domainmatrix import (
        DDM,
        DDMBadInputError, DDMDomainError, DDMShapeError,
        ddm_iadd, ddm_isub, ddm_ineg, ddm_imatmul,
        ddm_irref, ddm_idet, ddm_iinv, ddm_ilu, ddm_ilu_split,
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
    A = DDM([], (0, 0), ZZ)
    assert A.det() == ZZ(1)

    A = DDM([[ZZ(2)]], (1, 1), ZZ)
    assert A.det() == ZZ(2)

    A = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    assert A.det() == ZZ(-2)

    A = DDM([[ZZ(1), ZZ(2), ZZ(3)], [ZZ(1), ZZ(2), ZZ(4)], [ZZ(1), ZZ(2), ZZ(5)]], (3, 3), ZZ)
    assert A.det() == ZZ(0)

    A = DDM([[QQ(1, 2), QQ(1, 2)], [QQ(1, 3), QQ(1, 4)]], (2, 2), QQ)
    assert A.det() == QQ(-1, 24)

    A = DDM([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    raises(DDMShapeError, lambda: A.det())

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

    A = []
    Ar = []
    pivots = []
    assert ddm_irref(A) == pivots
    assert A == Ar

    A = [[QQ(0), QQ(1)], [QQ(1), QQ(1)]]
    Ar = [[QQ(1), QQ(0)], [QQ(0), QQ(1)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    A = [[QQ(1), QQ(2), QQ(1)], [QQ(3), QQ(4), QQ(1)]]
    Ar = [[QQ(1), QQ(0), QQ(-1)], [QQ(0), QQ(1), QQ(1)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    A = [[QQ(3), QQ(4), QQ(1)], [QQ(1), QQ(2), QQ(1)]]
    Ar = [[QQ(1), QQ(0), QQ(-1)], [QQ(0), QQ(1), QQ(1)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    A = [[QQ(1), QQ(0)], [QQ(1), QQ(3)], [QQ(0), QQ(1)]]
    Ar = [[QQ(1), QQ(0)], [QQ(0), QQ(1)], [QQ(0), QQ(0)]]
    pivots = [0, 1]
    assert ddm_irref(A) == pivots
    assert A == Ar

    A = [[QQ(1), QQ(0), QQ(1)], [QQ(3), QQ(0), QQ(1)]]
    Ar = [[QQ(1), QQ(0), QQ(0)], [QQ(0), QQ(0), QQ(1)]]
    pivots = [0, 2]
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
