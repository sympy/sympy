from sympy import PurePoly, Symbol
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
