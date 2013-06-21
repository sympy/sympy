from sympy.matrices.sparsetools import _doktocsr, _csrtodok, _mulspvec
from sympy import SparseMatrix, Matrix
from sympy import ZZ, QQ


def test_doktocsr():
    a = SparseMatrix([[1, 2, 0, 0], [0, 3, 9, 0], [0, 1, 4, 0]])
    b = SparseMatrix(4, 6, [10, 20, 0, 0, 0, 0, 0, 30, 0, 40, 0, 0, 0, 0, 50,
        60, 70, 0, 0, 0, 0, 0, 0, 80])
    c = SparseMatrix(4, 4, [0, 0, 0, 0, 0, 12, 0, 2, 15, 0, 12, 0, 0, 0, 0, 4])
    d = SparseMatrix(10, 10, {(1, 1): 12, (3, 5): 7, (7, 8): 12})
    e = SparseMatrix([[0, 0, 0], [1, 0, 2], [3, 0, 0]])
    f = SparseMatrix(7, 8, {(2, 3): 5, (4, 5):12})
    assert _doktocsr(a) == [[1, 2, 3, 9, 1, 4], [0, 1, 1, 2, 1, 2],
        [0, 2, 4, 6], [3, 4]]
    assert _doktocsr(b) == [[10, 20, 30, 40, 50, 60, 70, 80],
        [0, 1, 1, 3, 2, 3, 4, 5], [0, 2, 4, 7, 8], [4, 6]]
    assert _doktocsr(c) == [[12, 2, 15, 12, 4], [1, 3, 0, 2, 3],
        [0, 0, 2, 4, 5], [4, 4]]
    assert _doktocsr(d) == [[12, 7, 12], [1, 5, 8],
        [0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3], [10, 10]]
    assert _doktocsr(e) == [[1, 2, 3], [0, 2, 0], [0, 0, 2, 3], [3, 3]]
    assert _doktocsr(f) == [[5, 12], [3, 5], [0, 0, 0, 1, 1, 2, 2, 2], [7, 8]]


def test_csrtodok():
    h = [[5, 7, 5], [2, 1, 3], [0, 1, 1, 3], [3, 4]]
    g = [[12, 5, 4], [2, 4, 2], [0, 1, 2, 3], [3, 7]]
    i = [[1, 3, 12], [0, 2, 4], [0, 2, 3], [2, 5]]
    j = [[11, 15, 12, 15], [2, 4, 1, 2], [0, 1, 1, 2, 3, 4], [5, 8]]
    k = [[1, 3], [2, 1], [0, 1, 1, 2], [3, 3]]
    assert _csrtodok(h) == SparseMatrix(3, 4,
        {(0, 2): 5, (2, 1): 7, (2, 3): 5})
    assert _csrtodok(g) == SparseMatrix(3, 7,
        {(0, 2): 12, (1, 4): 5, (2, 2): 4})
    assert _csrtodok(i) == SparseMatrix([[1, 0, 3, 0, 0], [0, 0, 0, 0, 12]])
    assert _csrtodok(j) == SparseMatrix(5, 8,
        {(0, 2): 11, (2, 4): 15, (3, 1): 12, (4, 2): 15})
    assert _csrtodok(k) == SparseMatrix(3, 3, {(0, 2): 1, (2, 1): 3})


def test_mulspvec():
    i = [[ZZ(1), ZZ(3), ZZ(12)], [0, 2, 4], [0, 2, 3], [2, 5]]
    j = [[ZZ(1), ZZ(2), ZZ(3), ZZ(4), ZZ(5), ZZ(6), ZZ(7), ZZ(8)], [0, 1, 1, 3, 2, 3, 4, 5], [0, 2, 4, 7, 8], [4, 6]]
    k = [[QQ(12, 15), QQ(14, 15), QQ(2, 3), QQ(17, 3)], [0, 2, 1, 3], [0, 2, 3, 4], [3, 4]]
    m = Matrix([[ZZ(12)], [ZZ(25)], [ZZ(2)], [ZZ(7)], [ZZ(8)]])
    n = Matrix([[ZZ(4)], [ZZ(7)], [ZZ(12)], [ZZ(14)], [ZZ(17)], [ZZ(2)]])
    o = Matrix([[QQ(3,2)], [QQ(7, 2)], [QQ(12, 7)], [QQ(15, 16)]])
    assert _mulspvec(i, m, ZZ) == SparseMatrix(2, 1, {(0, 0): 18, (1, 0): 96})
    assert _mulspvec(j, n, ZZ) == SparseMatrix(4, 1, {(0, 0): 18, (1, 0): 77, (2, 0): 263, (3, 0): 16})
    #assert _mulspvec(k, o, QQ) == SparseMatrix(3, 1, {(0, 0)})
