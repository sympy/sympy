from sympy import S, PurePoly, Symbol
from sympy.matrices import NonSquareMatrixError, Matrix, NPMatrix
from sympy.utilities.pytest import raises, XFAIL, skip, warns_deprecated_sympy

try:
    import numpy as np

    skipcheck = lambda: None

except ImportError:
    skipcheck = lambda: skip('NumPy must be available to test NumPyMatrix')


x = Symbol('x')


def test_creation():
    skipcheck()

    a = np.array([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(a)
    assert M == Matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, a)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])

    m = np.matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(m)
    assert M == Matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, m)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])

    a = np.array([[1, 2], [3, 4]], dtype='i4')
    M = NPMatrix(m)
    assert M.dtype == np.int32
    M = NPMatrix(m, dtype='f4')
    assert M.dtype == np.float32


def test_as_wrapper():
    skipcheck()

    a = np.array([[1, 2], [3, 4]])
    M = NPMatrix(a, copy=False)
    assert M[0] == 1
    a[0, 0] = -1
    assert M[0] == -1

    m = np.matrix([[1, 2], [3, 4]])
    M = NPMatrix(m, copy=False)
    assert M[0] == 1
    m[0, 0] = -1
    assert M[0] == -1


def test_as_reshaping_wrapper():
    skipcheck()

    a = np.array([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, a, copy=False)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])
    a[0, 0] = -1
    assert M[0] == -1

    m = np.matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, m, copy=False)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])
    m[0, 0] = -1
    assert M[0] == -1


def test_adjugate():
    skipcheck()

    M = NPMatrix(4, 4, [0, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])
    assert M.adjugate() == Matrix([
        [   4,  -8,  4,   0],
        [  76, -68, -8,  24],
        [-122, 142,  4, -48],
        [  48, -72,  0,  24]])
    M = NPMatrix(8, 2, M)
    raises(NonSquareMatrixError, lambda: M.adjugate())


def test_charpoly():
    skipcheck()

    M = NPMatrix([[1, 2], [3, 4]])
    assert M.charpoly(x='x') == PurePoly(x**2 - 5*x - 2, x, domain='RR')
    M = NPMatrix(4, 1, M)
    raises(NonSquareMatrixError, lambda: M.charpoly())


def test_cofactor():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.cofactor(1, 1).epsilon_eq(-12)


def test_cofactor_matrix():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.cofactor_matrix() == Matrix([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])


def test_det():
    skipcheck()

    M = NPMatrix([[1, 2], [3, 4]])
    assert M.det().epsilon_eq(-2)
    assert M.det_bareiss().epsilon_eq(-2)
    assert M.det_berkowitz().epsilon_eq(-2)
    assert M.det_LU().epsilon_eq(-2)

    M = NPMatrix(4, 1, M)
    raises(NonSquareMatrixError, lambda: M.det())
    raises(NonSquareMatrixError, lambda: M.det_bareiss())
    raises(NonSquareMatrixError, lambda: M.det_berkowitz())
    raises(NonSquareMatrixError, lambda: M.det_LU())


def test_minor():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.cofactor(1, 1).epsilon_eq(-12)


def test_minor_submatrix():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.minor_submatrix(1, 1) == Matrix([[1, 3], [7, 9]])


def test_echelon():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.is_echelon is False
    assert M.echelon_form() == Matrix([[7, 8, 9], [0, 6, 12], [0, 0, 0]])


def test_rank():
    skipcheck()

    M = NPMatrix([[S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I'),S('3/64+13/64*I'),S('-23/32-59/256*I'),S('15/128-3/32*I'),S('19/256+551/1024*I')],
        [S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I'),S('-219/128+115/256*I'),S('6301/4096-6609/1024*I'),S('119/128+143/128*I'),S('-10879/2048+4343/4096*I'),S('129/256-549/512*I'),S('42533/16384+29103/8192*I')],
        [S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I'),S('3/64+13/64*I'),S('-23/32-59/256*I')],
        [S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I'),S('-219/128+115/256*I'),S('6301/4096-6609/1024*I'),S('119/128+143/128*I'),S('-10879/2048+4343/4096*I')],
        [S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I')],
        [S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I'),S('-219/128+115/256*I'),S('6301/4096-6609/1024*I')],
        [S('-4'),S('9-5*I'),S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I')],
        [S('-2*I'),S('119/8+29/4*I'),S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I')],
        [S('0'),S('-6'),S('-4'),S('9-5*I'),S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I')],
        [S('1'),S('-9/4+3*I'),S('-2*I'),S('119/8+29/4*I'),S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I')],
        [S('0'),S('-4*I'),S('0'),S('-6'),S('-4'),S('9-5*I'),S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I')],
        [S('0'),S('1/4+1/2*I'),S('1'),S('-9/4+3*I'),S('-2*I'),S('119/8+29/4*I'),S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I')]])
    assert M.rank() == 8


def test_rref():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.is_echelon is False
    assert M.rref() == (Matrix([[1, 0, -1], [0, 1, 2], [0, 0, 0]]), (0, 1))
