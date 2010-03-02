from sympy.mpmath import *

def test_matrix_basic():
    A1 = matrix(3)
    for i in xrange(3):
        A1[i,i] = 1
    assert A1 == eye(3)
    assert A1 == matrix(A1)
    A2 = matrix(3, 2)
    assert not A2._matrix__data
    A3 = matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert list(A3) == range(1, 10)
    A3[1,1] = 0
    assert not (1, 1) in A3._matrix__data
    A4 = matrix([[1, 2, 3], [4, 5, 6]])
    A5 = matrix([[6, -1], [3, 2], [0, -3]])
    assert A4 * A5 == matrix([[12, -6], [39, -12]])
    assert A1 * A3 == A3 * A1 == A3
    try:
        A2 * A2
        assert False
    except ValueError:
        pass
    l = [[10, 20, 30], [40, 0, 60], [70, 80, 90]]
    A6 = matrix(l)
    assert A6.tolist() == l
    assert A6 == eval(repr(A6))
    A6 = matrix(A6, force_type=float)
    assert A6 == eval(repr(A6))
    assert A6*1j == eval(repr(A6*1j))
    assert A3 * 10 == 10 * A3 == A6
    assert A2.rows == 3
    assert A2.cols == 2
    A3.rows = 2
    A3.cols = 2
    assert len(A3._matrix__data) == 3
    assert A4 + A4 == 2*A4
    try:
        A4 + A2
    except ValueError:
        pass
    assert sum(A1 - A1) == 0
    A7 = matrix([[1, 2], [3, 4], [5, 6], [7, 8]])
    x = matrix([10, -10])
    assert A7*x == matrix([-10, -10, -10, -10])
    A8 = ones(5)
    assert sum((A8 + 1) - (2 - zeros(5))) == 0
    assert (1 + ones(4)) / 2 - 1 == zeros(4)
    assert eye(3)**10 == eye(3)
    try:
        A7**2
        assert False
    except ValueError:
        pass
    A9 = randmatrix(3)
    A10 = matrix(A9)
    A9[0,0] = -100
    assert A9 != A10
    A11 = matrix(randmatrix(2, 3), force_type=mpi)
    for a in A11:
        assert isinstance(a, mpi)
    assert nstr(A9)

def test_matrix_power():
    A = matrix([[1, 2], [3, 4]])
    assert A**2 == A*A
    assert A**3 == A*A*A
    assert A**-1 == inverse(A)
    assert A**-2 == inverse(A*A)

def test_matrix_transform():
    A = matrix([[1, 2], [3, 4], [5, 6]])
    assert A.T == A.transpose() == matrix([[1, 3, 5], [2, 4, 6]])
    swap_row(A, 1, 2)
    assert A == matrix([[1, 2], [5, 6], [3, 4]])
    l = [1, 2]
    swap_row(l, 0, 1)
    assert l == [2, 1]
    assert extend(eye(3), [1,2,3]) == matrix([[1,0,0,1],[0,1,0,2],[0,0,1,3]])

def test_matrix_conjugate():
    A = matrix([[1 + j, 0], [2, j]])
    assert A.conjugate() == matrix([[mpc(1, -1), 0], [2, mpc(0, -1)]])
    assert A.transpose_conj() == A.H == matrix([[mpc(1, -1), 2],
                                                [0, mpc(0, -1)]])

def test_matrix_creation():
    assert diag([1, 2, 3]) == matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    A1 = ones(2, 3)
    assert A1.rows == 2 and A1.cols == 3
    for a in A1:
        assert a == 1
    A2 = zeros(3, 2)
    assert A2.rows == 3 and A2.cols == 2
    for a in A2:
        assert a == 0
    assert randmatrix(10) != randmatrix(10)
    one = mpf(1)
    assert hilbert(3) == matrix([[one, one/2, one/3],
                                 [one/2, one/3, one/4],
                                 [one/3, one/4, one/5]])

def test_norms():
    # matrix norms
    A = matrix([[1, -2], [-3, -1], [2, 1]])
    assert mnorm(A,1) == 6
    assert mnorm(A,inf) == 4
    assert mnorm(A,'F') == sqrt(20)
    # vector norms
    assert norm(-3) == 3
    x = [1, -2, 7, -12]
    assert norm(x, 1) == 22
    assert round(norm(x, 2), 10) == 14.0712472795
    assert round(norm(x, 10), 10) == 12.0054633727
    assert norm(x, inf) == 12

def test_vector():
    x = matrix([0, 1, 2, 3, 4])
    assert x == matrix([[0], [1], [2], [3], [4]])
    assert x[3] == 3
    assert len(x._matrix__data) == 4
    assert list(x) == range(5)
    x[0] = -10
    x[4] = 0
    assert x[0] == -10
    assert len(x) == len(x.T) == 5
    assert x.T*x == matrix([[114]])

def test_matrix_copy():
    A = ones(6)
    B = A.copy()
    assert A == B
    B[0,0] = 0
    assert A != B

def test_matrix_numpy():
    try:
        import numpy
    except ImportError:
        return
    l = [[1, 2], [3, 4], [5, 6]]
    a = numpy.matrix(l)
    assert matrix(l) == matrix(a)

