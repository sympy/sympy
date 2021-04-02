from sympy.testing.pytest import raises

from sympy.core.numbers import Rational
from sympy.functions import sqrt

from sympy.matrices.common import (NonInvertibleMatrixError,
    NonSquareMatrixError, ShapeError)
from sympy.matrices.dense import Matrix
from sympy.polys import ZZ, QQ

from sympy.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.matrices.exceptions import (DDMBadInputError, DDMDomainError,
        DDMShapeError, DDMFormatError)
from sympy.polys.matrices.ddm import DDM
from sympy.polys.matrices.sdm import SDM


def test_DomainMatrix_init():
    lol = [[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]]
    dod = {0: {0: ZZ(1), 1:ZZ(2)}, 1: {0:ZZ(3), 1:ZZ(4)}}
    ddm = DDM(lol, (2, 2), ZZ)
    sdm = SDM(dod, (2, 2), ZZ)

    A = DomainMatrix(lol, (2, 2), ZZ)
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    A = DomainMatrix(dod, (2, 2), ZZ)
    assert A.rep == sdm
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    raises(TypeError, lambda: DomainMatrix(ddm, (2, 2), ZZ))
    raises(TypeError, lambda: DomainMatrix(sdm, (2, 2), ZZ))
    raises(TypeError, lambda: DomainMatrix(Matrix([[1]]), (1, 1), ZZ))

    for fmt, rep in [('sparse', sdm), ('dense', ddm)]:
        A = DomainMatrix(lol, (2, 2), ZZ, fmt=fmt)
        assert A.rep == rep
        A = DomainMatrix(dod, (2, 2), ZZ, fmt=fmt)
        assert A.rep == rep

    raises(ValueError, lambda: DomainMatrix(lol, (2, 2), ZZ, fmt='invalid'))

    raises(DDMBadInputError, lambda: DomainMatrix([[ZZ(1), ZZ(2)]], (2, 2), ZZ))


def test_DomainMatrix_from_rep():
    ddm = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A = DomainMatrix.from_rep(ddm)
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    sdm = SDM({0: {0: ZZ(1), 1:ZZ(2)}, 1: {0:ZZ(3), 1:ZZ(4)}}, (2, 2), ZZ)
    A = DomainMatrix.from_rep(sdm)
    assert A.rep == sdm
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    A = DomainMatrix([[ZZ(1)]], (1, 1), ZZ)
    raises(TypeError, lambda: DomainMatrix.from_rep(A))


def test_DomainMatrix_from_list_sympy():
    ddm = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A = DomainMatrix.from_list_sympy(2, 2, [[1, 2], [3, 4]])
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    K = QQ.algebraic_field(sqrt(2))
    ddm = DDM(
        [[K.convert(1 + sqrt(2)), K.convert(2 + sqrt(2))],
         [K.convert(3 + sqrt(2)), K.convert(4 + sqrt(2))]],
        (2, 2),
        K
    )
    A = DomainMatrix.from_list_sympy(
        2, 2, [[1 + sqrt(2), 2 + sqrt(2)], [3 + sqrt(2), 4 + sqrt(2)]],
        extension=True)
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == K


def test_DomainMatrix_from_Matrix():
    ddm = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A = DomainMatrix.from_Matrix(Matrix([[1, 2], [3, 4]]))
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == ZZ

    K = QQ.algebraic_field(sqrt(2))
    ddm = DDM(
        [[K.convert(1 + sqrt(2)), K.convert(2 + sqrt(2))],
         [K.convert(3 + sqrt(2)), K.convert(4 + sqrt(2))]],
        (2, 2),
        K
    )
    A = DomainMatrix.from_Matrix(
        Matrix([[1 + sqrt(2), 2 + sqrt(2)], [3 + sqrt(2), 4 + sqrt(2)]]),
        extension=True)
    assert A.rep == ddm
    assert A.shape == (2, 2)
    assert A.domain == K


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


def test_DomainMatrix_to_sparse():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A_sparse = A.to_sparse()
    assert A_sparse.rep == {0: {0: 1, 1: 2}, 1: {0: 3, 1: 4}}


def test_DomainMatrix_to_dense():
    A = DomainMatrix({0: {0: 1, 1: 2}, 1: {0: 3, 1: 4}}, (2, 2), ZZ)
    A_dense = A.to_dense()
    assert A_dense.rep == DDM([[1, 2], [3, 4]], (2, 2), ZZ)


def test_DomainMatrix_unify():
    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    assert Az.unify(Az) == (Az, Az)
    assert Az.unify(Aq) == (Aq, Aq)
    assert Aq.unify(Az) == (Aq, Aq)
    assert Aq.unify(Aq) == (Aq, Aq)

    As = DomainMatrix({0: {1: ZZ(1)}, 1:{0:ZZ(2)}}, (2, 2), ZZ)
    Ad = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)

    assert As.unify(As) == (As, As)
    assert Ad.unify(Ad) == (Ad, Ad)

    Bs, Bd = As.unify(Ad, fmt='dense')
    assert Bs.rep == DDM([[0, 1], [2, 0]], (2, 2), ZZ)
    assert Bd.rep == DDM([[1, 2],[3, 4]], (2, 2), ZZ)

    Bs, Bd = As.unify(Ad, fmt='sparse')
    assert Bs.rep == SDM({0: {1: 1}, 1: {0: 2}}, (2, 2), ZZ)
    assert Bd.rep == SDM({0: {0: 1, 1: 2}, 1: {0: 3, 1: 4}}, (2, 2), ZZ)

    raises(ValueError, lambda: As.unify(Ad, fmt='invalid'))


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
    raises(DDMShapeError, lambda: A1 + A2)
    raises(DDMShapeError, lambda: A2 + A1)
    raises(DDMShapeError, lambda: A1.add(A2))
    raises(DDMShapeError, lambda: A2.add(A1))

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    Asum = DomainMatrix([[QQ(2), QQ(4)], [QQ(6), QQ(8)]], (2, 2), QQ)
    assert Az + Aq == Asum
    assert Aq + Az == Asum
    raises(DDMDomainError, lambda: Az.add(Aq))
    raises(DDMDomainError, lambda: Aq.add(Az))

    As = DomainMatrix({0: {1: ZZ(1)}, 1: {0: ZZ(2)}}, (2, 2), ZZ)
    Ad = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)

    Asd = As + Ad
    Ads = Ad + As
    assert Asd == DomainMatrix([[1, 3], [5, 4]], (2, 2), ZZ)
    assert Asd.rep == DDM([[1, 3], [5, 4]], (2, 2), ZZ)
    assert Ads == DomainMatrix([[1, 3], [5, 4]], (2, 2), ZZ)
    assert Ads.rep == DDM([[1, 3], [5, 4]], (2, 2), ZZ)
    raises(DDMFormatError, lambda: As.add(Ad))


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
    raises(DDMShapeError, lambda: A1 - A2)
    raises(DDMShapeError, lambda: A2 - A1)
    raises(DDMShapeError, lambda: A1.sub(A2))
    raises(DDMShapeError, lambda: A2.sub(A1))

    Az = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    Aq = DomainMatrix([[QQ(1), QQ(2)], [QQ(3), QQ(4)]], (2, 2), QQ)
    Adiff = DomainMatrix([[QQ(0), QQ(0)], [QQ(0), QQ(0)]], (2, 2), QQ)
    assert Az - Aq == Adiff
    assert Aq - Az == Adiff
    raises(DDMDomainError, lambda: Az.sub(Aq))
    raises(DDMDomainError, lambda: Aq.sub(Az))

    As = DomainMatrix({0: {1: ZZ(1)}, 1: {0: ZZ(2)}}, (2, 2), ZZ)
    Ad = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)

    Asd = As - Ad
    Ads = Ad - As
    assert Asd == DomainMatrix([[-1, -1], [-1, -4]], (2, 2), ZZ)
    assert Asd.rep == DDM([[-1, -1], [-1, -4]], (2, 2), ZZ)
    assert Asd == -Ads
    assert Asd.rep == -Ads.rep


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
    Aprod = DomainMatrix([[QQ(7), QQ(10)], [QQ(15), QQ(22)]], (2, 2), QQ)
    assert Az * Aq == Aprod
    assert Aq * Az == Aprod
    raises(DDMDomainError, lambda: Az.matmul(Aq))
    raises(DDMDomainError, lambda: Aq.matmul(Az))

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    AA = DomainMatrix([[ZZ(2), ZZ(4)], [ZZ(6), ZZ(8)]], (2, 2), ZZ)
    x = ZZ(2)
    assert A * x == x * A == A.mul(x) == AA

    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    AA = DomainMatrix([[ZZ(0), ZZ(0)], [ZZ(0), ZZ(0)]], (2, 2), ZZ)
    x = ZZ(0)
    assert A * x == x * A == A.mul(x) == AA

    As = DomainMatrix({0: {1: ZZ(1)}, 1: {0: ZZ(2)}}, (2, 2), ZZ)
    Ad = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)

    Asd = As * Ad
    Ads = Ad * As
    assert Asd == DomainMatrix([[3, 4], [2, 4]], (2, 2), ZZ)
    assert Asd.rep == DDM([[3, 4], [2, 4]], (2, 2), ZZ)
    assert Ads == DomainMatrix([[4, 1], [8, 3]], (2, 2), ZZ)
    assert Ads.rep == DDM([[4, 1], [8, 3]], (2, 2), ZZ)


def test_DomainMatrix_pow():
    eye = DomainMatrix.eye(2, ZZ)
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    A2 = DomainMatrix([[ZZ(7), ZZ(10)], [ZZ(15), ZZ(22)]], (2, 2), ZZ)
    A3 = DomainMatrix([[ZZ(37), ZZ(54)], [ZZ(81), ZZ(118)]], (2, 2), ZZ)
    assert A**0 == A.pow(0) == eye
    assert A**1 == A.pow(1) == A
    assert A**2 == A.pow(2) == A2
    assert A**3 == A.pow(3) == A3

    raises(TypeError, lambda: A ** Rational(1, 2))
    raises(NotImplementedError, lambda: A ** -1)
    raises(NotImplementedError, lambda: A.pow(-1))

    A = DomainMatrix.zeros((2, 1), ZZ)
    raises(NonSquareMatrixError, lambda: A ** 1)


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


def test_DomainMatrix_nullspace():
    A = DomainMatrix([[QQ(1), QQ(1)], [QQ(1), QQ(1)]], (2, 2), ZZ)
    Anull = DomainMatrix([[QQ(-1), QQ(1)]], (1, 2), ZZ)
    assert A.nullspace() == Anull


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


def test_DomainMatrix_eye():
    A = DomainMatrix.eye(3, QQ)
    assert A.rep == SDM.eye(3, QQ)
    assert A.shape == (3, 3)
    assert A.domain == QQ


def test_DomainMatrix_zeros():
    A = DomainMatrix.zeros((1, 2), QQ)
    assert A.rep == SDM.zeros((1, 2), QQ)
    assert A.shape == (1, 2)
    assert A.domain == QQ


def test_DomainMatrix_hstack():
    A = DomainMatrix([[ZZ(1)], [ZZ(2)]], (2, 1), ZZ)
    B = DomainMatrix([[QQ(3), QQ(4)], [QQ(5), QQ(6)]], (2, 2), QQ)
    AB = DomainMatrix([[QQ(1), QQ(3), QQ(4)], [QQ(2), QQ(5), QQ(6)]], (2, 3), QQ)
    assert A.hstack(B) == AB
