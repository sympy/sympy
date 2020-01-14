from sympy import (Float, I, N, PurePoly, S, Symbol, oo, nan, sqrt, sympify,
    Matrix, ImmutableDenseMatrix, NumPyMatrix, NonSquareMatrixError,
    diag, eye, zeros)
from sympy.matrices.matrices import MatrixError
from sympy.testing.pytest import raises, skip

try:
    import numpy as np

    skipcheck = lambda: None

except ImportError:
    skipcheck = lambda: skip('NumPy must be available to test NumPyMatrix')

x = Symbol('x')


def _fcmp(a, b):
    return Float(N(a)).epsilon_eq(b, epsilon=1e-12)

def _fcmpseq(seq1, seq2):
    return all(_fcmp(a, b) for a, b in zip(seq1, seq2))

def _fcmpseqsort(seq1, seq2):
    return all(_fcmp(a, b) for a, b in zip(sorted(seq1), sorted(seq2)))


# creation

def test_creation():
    skipcheck()

    a = np.array([[1, 2, 3], [4, 5, 6]])
    M = NumPyMatrix(a)
    assert M == Matrix([[1, 2, 3], [4, 5, 6]])
    M = NumPyMatrix(3, 2, a)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])

    a = np.array([[1, 2], [3, 4]], dtype='i4')
    M = NumPyMatrix(a)
    assert M.dtype == np.int32
    M = NumPyMatrix(a, dtype='f4')
    assert M.dtype == np.float32

    m = np.matrix([[1, 2, 3], [4, 5, 6]])
    M = NumPyMatrix(m)
    assert M == Matrix([[1, 2, 3], [4, 5, 6]])
    M = NumPyMatrix(3, 2, m)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])

    m = np.matrix([[1, 2, 3], [4, 5, 6]], dtype='i4')
    M = NumPyMatrix(m)
    assert M.dtype == np.int32
    M = NumPyMatrix(m, dtype='f4')
    assert M.dtype == np.float32

def test_detect_dtype():
    skipcheck()

    assert NumPyMatrix([[1]]).dtype == np.int64
    assert NumPyMatrix([[1.0]]).dtype == np.float64
    assert NumPyMatrix([[oo]]).dtype == np.float64
    assert NumPyMatrix([[nan]]).dtype == np.float64
    assert NumPyMatrix([[I]]).dtype == np.complex128

def test_as_wrapper():
    skipcheck()

    a = np.array([[1, 2], [3, 4]])
    M = NumPyMatrix(a, copy=False)
    assert M[0] == 1
    a[0, 0] = -1
    assert M[0] == -1

    m = np.matrix([[1, 2], [3, 4]])
    M = NumPyMatrix(m, copy=False)
    assert M[0] == 1
    m[0, 0] = -1
    assert M[0] == -1

def test_as_reshaping_wrapper():
    skipcheck()

    a = np.array([[1, 2, 3], [4, 5, 6]])
    M = NumPyMatrix(3, 2, a, copy=False)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])
    a[0, 0] = -1
    assert M[0] == -1

    m = np.matrix([[1, 2, 3], [4, 5, 6]])
    M = NumPyMatrix(3, 2, m, copy=False)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])
    m[0, 0] = -1
    assert M[0] == -1

def test_sympify():
    skipcheck()

    a = np.array([[1, 2], [3, 4]])
    assert isinstance(sympify(NumPyMatrix(a)), ImmutableDenseMatrix)


# basic operations

def test_array():
    skipcheck()

    a = np.array([[1, 2], [3, 4]])
    M = NumPyMatrix(a, copy=False)
    assert M.__array__() is a
    assert M.__array__('f4') is not a

def test_abs():
    skipcheck()

    M = NumPyMatrix([[1, 1j], [-1j, 1+1j]])
    A = abs(M)

    assert A.dtype == np.float64
    assert _fcmpseq(A, (1, 1, 1, sqrt(2)))

def test_neg():
    skipcheck()

    M = NumPyMatrix([[1, -1], [-1, 2]], dtype='i8')
    A = -M

    assert A.dtype == np.int64
    assert _fcmpseq(A, (-1, 1, 1, -2))

def test_add():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]], dtype='i8')
    A = Matrix([[2, 4], [6, 8]])

    R = M + M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M + M._arr
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M._arr + M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M + NumPyMatrix(M, dtype='f8')
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.float64
    assert _fcmpseq(R, A)

    R = M + NumPyMatrix(M, dtype='c16')
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.complex128
    assert _fcmpseq(R, A)

    raises(TypeError, lambda: M + 1)
    raises(TypeError, lambda: 1 + M)

def test_sub():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]], dtype='i8')
    A = zeros(2)

    R = M - M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M - M._arr
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M._arr - M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

def test_mul():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]], dtype='i8')
    A = Matrix([[7, 10], [15, 22]])

    R = M * M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M * M._arr
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M._arr * M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M @ M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M @ M._arr
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M._arr @ M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M * NumPyMatrix(M, dtype='f8')
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.float64
    assert _fcmpseq(R, A)

    R = M * NumPyMatrix(M, dtype='c16')
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.complex128
    assert _fcmpseq(R, A)

    A = Matrix([[2, 4], [6, 8]])

    R = M * 2
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = 2 * M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert R == A

    R = M * 2.0
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.float64
    assert _fcmpseq(R, A)

    R = 2.0 * M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.float64
    assert _fcmpseq(R, A)

    R = M * (2+0j)
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.complex128
    assert _fcmpseq(R, A)

    R = (2+0j) * M
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.complex128
    assert _fcmpseq(R, A)

    A = Matrix([[x, 2*x], [3*x, 4*x]])

    R = x * M
    assert isinstance(R, Matrix)
    assert R == A

    R = M * x
    assert isinstance(R, Matrix)
    assert R == A

def test_div():
    skipcheck()

    M = NumPyMatrix([[2, 4], [6, 8]], dtype='i8')
    A = Matrix([[1, 2], [3, 4]])

    R = M / 2
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.float64
    assert _fcmpseq(R, A)

    R = M / 2.0
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.float64
    assert _fcmpseq(R, A)

    R = M / (2+0j)
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.complex128
    assert _fcmpseq(R, A)

    R = M / x
    assert isinstance(R, Matrix)
    assert R == Matrix([[2/x, 4/x], [6/x, 8/x]])

def test_pow():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]], dtype='i8')
    A = Matrix([[7, 10], [15, 22]])

    R = M**2
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert _fcmpseq(R, A)

    R = M**2.0
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert _fcmpseq(R, A)

    R = M**(2+0j)
    assert isinstance(R, NumPyMatrix)
    assert R.dtype == np.int64
    assert _fcmpseq(R, A)

    R = M**x
    assert isinstance(R, Matrix)
    assert R == Matrix(S('''[
        [-2*(5/2 - sqrt(33)/2)**x/((-3/2 + sqrt(33)/2)*(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2))) + 2*(5/2 + sqrt(33)/2)**x/((-sqrt(33)/2 - 3/2)*(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2))), -4*(5/2 - sqrt(33)/2)**x/((-3/2 + sqrt(33)/2)*(-sqrt(33)/2 - 3/2)*(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2))) + 4*(5/2 + sqrt(33)/2)**x/((-3/2 + sqrt(33)/2)*(-sqrt(33)/2 - 3/2)*(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2)))],
        [                                                 (5/2 - sqrt(33)/2)**x/(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2)) - (5/2 + sqrt(33)/2)**x/(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2)),                                          2*(5/2 - sqrt(33)/2)**x/((-sqrt(33)/2 - 3/2)*(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2))) - 2*(5/2 + sqrt(33)/2)**x/((-3/2 + sqrt(33)/2)*(-2/(-3/2 + sqrt(33)/2) + 2/(-sqrt(33)/2 - 3/2)))]]'''))

def test_coercion():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]], dtype='i8')
    assert isinstance(M * x, Matrix)
    assert isinstance(x * M, Matrix)
    assert isinstance(M / x, Matrix)
    assert isinstance(M + Matrix([[x, x], [x, x]]), Matrix)


# higher level functions

def test_adjugate():
    skipcheck()

    M = NumPyMatrix(4, 4, [0, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])
    A = M.adjugate()

    assert isinstance(A, NumPyMatrix)
    assert _fcmpseq(A, (4, -8, 4, 0, 76, -68, -8, 24, -122, 142, 4, -48, 48, -72, 0, 24))

    M = NumPyMatrix(8, 2, M)
    raises(NonSquareMatrixError, lambda: M.adjugate())

def test_charpoly():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]])
    assert M.charpoly(x='x') == PurePoly(x**2 - 5*x - 2, x, domain='RR')

    M = NumPyMatrix(4, 1, M)
    raises(NonSquareMatrixError, lambda: M.charpoly())

def test_cofactor():
    skipcheck()

    M = NumPyMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert _fcmp(M.cofactor(1, 1), -12)

def test_cofactor_matrix():
    skipcheck()

    M = NumPyMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = M.cofactor_matrix()

    assert isinstance(A, NumPyMatrix)
    assert _fcmpseq(A, (-3, 6, -3, 6, -12, 6, -3, 6, -3))

def test_det():
    skipcheck()

    M = NumPyMatrix([[1, 2], [3, 4]])
    assert _fcmp(M.det(), -2)
    assert _fcmp(M.det_bareiss(), -2)
    assert _fcmp(M.det_berkowitz(), -2)
    assert _fcmp(M.det_LU(), -2)

    M = NumPyMatrix(4, 1, M)
    raises(NonSquareMatrixError, lambda: M.det())
    raises(NonSquareMatrixError, lambda: M.det_bareiss())
    raises(NonSquareMatrixError, lambda: M.det_berkowitz())
    raises(NonSquareMatrixError, lambda: M.det_LU())

def test_minor():
    skipcheck()

    M = NumPyMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert _fcmp(M.cofactor(1, 1), -12)

def test_minor_submatrix():
    skipcheck()

    M = NumPyMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = M.minor_submatrix(1, 1)

    assert isinstance(A, NumPyMatrix)
    assert _fcmpseq(A, Matrix([[1, 3], [7, 9]]))

def test_echelon():
    skipcheck()

    M = NumPyMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype='i4') # float gives different result
    A = M.echelon_form()

    assert M.is_echelon is False
    assert isinstance(A, NumPyMatrix)
    assert _fcmpseq(A, Matrix([[1, 2, 3], [0, -3, -6], [0, 0, 0]]))

def test_rank():
    skipcheck()

    M = NumPyMatrix([[S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I'),S('3/64+13/64*I'),S('-23/32-59/256*I'),S('15/128-3/32*I'),S('19/256+551/1024*I')],
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

    M = NumPyMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = M.rref()

    assert isinstance(A[0], NumPyMatrix)
    assert _fcmpseq(A[0], Matrix([[1, 0, -1], [0, 1, 2], [0, 0, 0]]))
    assert A[1] == (0, 1)

def test_columnspace():
    skipcheck()

    M = NumPyMatrix([
        [ 1,  2,  0,  2,  5],
        [-2, -5,  1, -1, -8],
        [ 0, -3,  3,  4,  1],
        [ 3,  6,  0, -7,  2]])
    A = M.columnspace()

    assert all(isinstance(m, NumPyMatrix) for m in A)
    assert _fcmpseq(A[0], Matrix([1, -2, 0, 3]))
    assert _fcmpseq(A[1], Matrix([2, -5, -3, 6]))
    assert _fcmpseq(A[2], Matrix([2, -1, 4, -7]))

def test_rowspace():
    skipcheck()

    M = NumPyMatrix([
        [ 1,  2,  0,  2,  5],
        [-2, -5,  1, -1, -8],
        [ 0, -3,  3,  4,  1],
        [ 3,  6,  0, -7,  2]], dtype='i4') # float gives different result
    A = M.rowspace()

    assert all(isinstance(m, NumPyMatrix) for m in A)
    assert _fcmpseq(A[0], Matrix([[1, 2, 0, 2, 5]]))
    assert _fcmpseq(A[1], Matrix([[0, -1, 1, 3, 2]]))
    assert _fcmpseq(A[2], Matrix([[0, 0, 0, 5, 5]]))

def test_nullspace():
    skipcheck()

    M = NumPyMatrix([
        [ 1,  3, 0,  2,  6, 3, 1],
        [-2, -6, 0, -2, -8, 3, 1],
        [ 3,  9, 0,  0,  6, 6, 2],
        [-1, -3, 0,  1,  0, 9, 3]])
    A = M.nullspace()

    assert all(isinstance(m, NumPyMatrix) for m in A)
    assert _fcmpseq(A[0], Matrix([-3, 1, 0, 0, 0, 0, 0]))
    assert _fcmpseq(A[1], Matrix([0, 0, 1, 0, 0, 0, 0]))
    assert _fcmpseq(A[2], Matrix([-2, 0, 0, -2, 1, 0, 0]))
    assert _fcmpseq(A[3], Matrix([0, 0, 0, 0, 0, -1/3, 1]))

def test_orthogonalize():
    skipcheck()

    v = [Matrix([1, I]), Matrix([1, -I])]
    A = NumPyMatrix.orthogonalize(*v)

    assert all(isinstance(m, NumPyMatrix) for m in A)
    assert _fcmp(A[0][0], 1)
    assert _fcmp(A[0][1] / I, 1)
    assert _fcmp(A[1][0], 1)
    assert _fcmp(A[1][1] / I, -1)

def test_eigenvals():
    skipcheck()

    M = NumPyMatrix(3, 3, [0, 1, 1, 1, 0, 0, 1, 1, 1])
    assert M.eigenvals() == {2: 1, -1: 1, 0: 1}

    M = NumPyMatrix(3, 3, [1.5, 0.75, 1.75, 0, -0.75, 0.5, 0, 0, -0.25])
    e = M.eigenvals()
    assert _fcmpseqsort(e.keys(), (1.5, -0.75, -0.25))
    assert all(v == 1 for v in e.values())

    M = NumPyMatrix(3, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1])
    assert M.eigenvals() == {1: 3}
    assert M.eigenvals(multiple=True) == [1, 1, 1]

def test_eigenvects():
    skipcheck()

    M = Matrix(3, 3, [0, 1, 1, 1, 0, 0, 1, 1, 1])
    A = M.eigenvects()
    B = NumPyMatrix(M).eigenvects()

    assert _fcmp(B[0][0], A[0][0])
    assert _fcmp(B[1][0], A[1][0])
    assert _fcmp(B[2][0], A[2][0])
    assert B[0][1] == A[0][1]
    assert B[1][1] == A[1][1]
    assert B[2][1] == A[2][1]
    assert _fcmpseq(A[0][2][0], B[0][2][0])
    assert _fcmpseq(A[1][2][0], B[1][2][0])
    assert _fcmpseq(A[2][2][0], B[2][2][0])

    A = M.left_eigenvects()
    B = NumPyMatrix(M).left_eigenvects()

    assert _fcmp(B[0][0], A[0][0])
    assert _fcmp(B[1][0], A[1][0])
    assert _fcmp(B[2][0], A[2][0])
    assert B[0][1] == A[0][1]
    assert B[1][1] == A[1][1]
    assert B[2][1] == A[2][1]
    assert _fcmpseq(A[0][2][0], B[0][2][0])
    assert _fcmpseq(A[1][2][0], B[1][2][0])
    assert _fcmpseq(A[2][2][0], B[2][2][0])

def test_diagonalization():
    skipcheck ()

    M = NumPyMatrix([[1, 2+I], [2-I, 3]])
    assert M.is_diagonalizable()

    M = NumPyMatrix(3, 2, [-3, 1, -3, 20, 3, 10])
    assert not M.is_diagonalizable()
    assert not M.is_symmetric()
    raises(NonSquareMatrixError, lambda: M.diagonalize())

    # diagonalizable
    M = NumPyMatrix(diag(1, 2, 3))
    (P, D) = M.diagonalize()
    assert _fcmpseq(P, eye(3))
    assert _fcmpseq(D, M)

    M = NumPyMatrix(2, 2, [0, 1, 1, 0])
    assert M.is_symmetric()
    assert M.is_diagonalizable()
    (P, D) = M.diagonalize()
    assert _fcmpseq(P.inv() * M * P, D)

    M = NumPyMatrix(2, 2, [1, 0, 0, 3])
    assert M.is_symmetric()
    assert M.is_diagonalizable()
    (P, D) = M.diagonalize()
    assert _fcmpseq(P.inv() * M * P, D)
    assert _fcmpseq(P, eye(2))
    assert _fcmpseq(D, M)

    M = NumPyMatrix(2, 2, [1, 1, 0, 0])
    assert M.is_diagonalizable()
    (P, D) = M.diagonalize()
    assert _fcmpseq(P.inv() * M * P, D)

    M = NumPyMatrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
    assert M.is_diagonalizable()
    (P, D) = M.diagonalize()
    assert _fcmpseq(P.inv() * M * P, D)
    for i in P:
        assert i.as_numer_denom()[1] == 1

    M = NumPyMatrix(2, 2, [1, 0, 0, 0])
    assert M.is_diagonal()
    assert M.is_diagonalizable()
    (P, D) = M.diagonalize()
    assert _fcmpseq(P.inv() * M * P, D)
    assert _fcmpseq(P, Matrix([[0, 1], [1, 0]]))

    # diagonalizable, complex only
    M = NumPyMatrix(2, 2, [0, 1, -1, 0])
    assert not M.is_diagonalizable(True)
    raises(MatrixError, lambda: M.diagonalize(True))
    assert M.is_diagonalizable()
    (P, D) = M.diagonalize()
    assert P.inv() * M * P == D # _fcmpseq(P.inv() * M * P, D)

    # not diagonalizable
    M = NumPyMatrix(2, 2, [0, 1, 0, 0])
    assert not M.is_diagonalizable()
    raises(MatrixError, lambda: M.diagonalize())

    M = NumPyMatrix(3, 3, [-3, 1, -3, 20, 3, 10, 2, -2, 4])
    assert not M.is_diagonalizable()
    raises(MatrixError, lambda: M.diagonalize())

def test_definite():
    skipcheck ()

    # Examples from Gilbert Strang, "Introduction to Linear Algebra"
    # Positive definite matrices
    M = NumPyMatrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
    assert M.is_positive_definite == True
    assert M.is_positive_semidefinite == True
    assert M.is_negative_definite == False
    assert M.is_negative_semidefinite == False
    assert M.is_indefinite == False

    M = NumPyMatrix([[5, 4], [4, 5]])
    assert M.is_positive_definite == True
    assert M.is_positive_semidefinite == True
    assert M.is_negative_definite == False
    assert M.is_negative_semidefinite == False
    assert M.is_indefinite == False

    # Positive semidefinite matrices
    M = NumPyMatrix([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]])
    assert M.is_positive_definite == False
    assert M.is_positive_semidefinite == True
    assert M.is_negative_definite == False
    assert M.is_negative_semidefinite == False
    assert M.is_indefinite == False

    M = NumPyMatrix([[1, 2], [2, 4]])
    assert M.is_positive_definite == False
    assert M.is_positive_semidefinite == True
    assert M.is_negative_definite == False
    assert M.is_negative_semidefinite == False
    assert M.is_indefinite == False

    # Examples from Mathematica documentation
    # Non-hermitian positive definite matrices
    M = NumPyMatrix([[2, 3], [4, 8]])
    assert M.is_positive_definite == True
    assert M.is_positive_semidefinite == True
    assert M.is_negative_definite == False
    assert M.is_negative_semidefinite == False
    assert M.is_indefinite == False

    M = NumPyMatrix([[1, 2*I], [-I, 4]])
    assert M.is_positive_definite == True
    assert M.is_positive_semidefinite == True
    assert M.is_negative_definite == False
    assert M.is_negative_semidefinite == False
    assert M.is_indefinite == False

def test_jordan_form():
    skipcheck ()

    m = NumPyMatrix(3, 2, [-3, 1, -3, 20, 3, 10])
    raises(NonSquareMatrixError, lambda: m.jordan_form())

    # diagonalizable
    M = NumPyMatrix(3, 3, [7, -12, 6, 10, -19, 10, 12, -24, 13])
    N = NumPyMatrix(3, 3, [-1, 0, 0, 0, 1, 0, 0, 0, 1])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)
    assert _fcmpseq(N, M.diagonalize()[1])

    # Jordan cells
    # complexity: one of eigenvalues is zero
    M = NumPyMatrix(3, 3, [0, 1, 0, -4, 4, 0, -2, 1, 2])
    N = NumPyMatrix(3, 3, [2, 1, 0, 0, 2, 0, 0, 0, 2])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)

    # complexity: all of eigenvalues are equal
    M = NumPyMatrix(3, 3, [2, 6, -15, 1, 1, -5, 1, 2, -6])
    N = NumPyMatrix(3, 3, [-1, 1, 0, 0, -1, 0, 0, 0, -1])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)

    # complexity: two of eigenvalues are zero
    M = NumPyMatrix(3, 3, [4, -5, 2, 5, -7, 3, 6, -9, 4])
    N = NumPyMatrix(3, 3, [0, 1, 0, 0, 0, 0, 0, 0, 1])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)

    M = NumPyMatrix(4, 4, [6, 2, -8, -6, -3, 2, 9, 6, 2, -2, -8, -6, -1, 0, 3, 4])
    N = NumPyMatrix(4, 4, [-2, 0, 0, 0,
                         0, 2, 1, 0,
                         0, 0, 2, 0,
                         0, 0, 0, 2])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)

    M = NumPyMatrix(4, 4, [6, 5, -2, -3, -3, -1, 3, 3, 2, 1, -2, -3, -1, 1, 5, 5])
    N = NumPyMatrix(4, 4, [2, 1, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)

    M = NumPyMatrix(4, 4, [5, 4, 2, 1, 0, 1, -1, -1, -1, -1, 3, 0, 1, 1, -1, 2])
    assert not M.is_diagonalizable()
    N = NumPyMatrix(4, 4, [1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 1, 0, 0, 0, 4])
    _, J = M.jordan_form()
    assert _fcmpseq(N, J)

def test_singular_values():
    skipcheck ()

    from math import sqrt

    M = NumPyMatrix([[0, 1*I], [2, 0]])
    assert _fcmpseqsort(M.singular_values(), (2, 1))

    M = NumPyMatrix([[2, 4], [1, 3], [0, 0], [0, 0]])
    assert _fcmpseqsort(M.singular_values(), [sqrt(sqrt(221) + 15), sqrt(15 - sqrt(221))])
    assert _fcmpseqsort(M.T.singular_values(), [sqrt(sqrt(221) + 15), sqrt(15 - sqrt(221)), 0, 0])
