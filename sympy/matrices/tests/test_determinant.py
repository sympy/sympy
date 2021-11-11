import random
from sympy.core.numbers import I
from sympy.core.numbers import Rational
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys.polytools import Poly
from sympy.matrices import Matrix, eye, ones
from sympy.abc import x, y, z
from sympy.testing.pytest import raises
from sympy.matrices.common import NonSquareMatrixError
from sympy.functions.combinatorial.factorials import factorial, subfactorial


def test_determinant():
    M = Matrix()

    assert M.det() == 1
    # Evaluating these directly because they are never reached via M.det()
    assert M._eval_det_bareiss() == 1
    assert M._eval_det_berkowitz() == 1
    assert M._eval_det_lu() == 1

    M = Matrix([ [0] ])

    assert M.det() == 0
    assert M._eval_det_bareiss() == 0
    assert M._eval_det_berkowitz() == 0
    assert M._eval_det_lu() == 0

    M = Matrix([ [5] ])

    assert M.det() == 5
    assert M._eval_det_bareiss() == 5
    assert M._eval_det_berkowitz() == 5
    assert M._eval_det_lu() == 5

    M = Matrix(( (-3,  2),
                 ( 8, -5) ))

    assert M.det(method="domain-ge") == -1
    assert M.det(method="bareiss") == -1
    assert M.det(method="berkowitz") == -1
    assert M.det(method="lu") == -1

    M = Matrix(( (x,   1),
                 (y, 2*y) ))

    assert M.det(method="domain-ge") == 2*x*y - y
    assert M.det(method="bareiss") == 2*x*y - y
    assert M.det(method="berkowitz") == 2*x*y - y
    assert M.det(method="lu") == 2*x*y - y

    M = Matrix(( (1, 1, 1),
                 (1, 2, 3),
                 (1, 3, 6) ))

    assert M.det(method="domain-ge") == 1
    assert M.det(method="bareiss") == 1
    assert M.det(method="berkowitz") == 1
    assert M.det(method="lu") == 1

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="domain-ge") == -289
    assert M.det(method="bareiss") == -289
    assert M.det(method="berkowitz") == -289
    assert M.det(method="lu") == -289

    M = Matrix(( ( 1,  2,  3,  4),
                 ( 5,  6,  7,  8),
                 ( 9, 10, 11, 12),
                 (13, 14, 15, 16) ))

    assert M.det(method="domain-ge") == 0
    assert M.det(method="bareiss") == 0
    assert M.det(method="berkowitz") == 0
    assert M.det(method="lu") == 0

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="domain-ge") == 275
    assert M.det(method="bareiss") == 275
    assert M.det(method="berkowitz") == 275
    assert M.det(method="lu") == 275

    M = Matrix(( ( 3,  0,  0, 0),
                 (-2,  1,  0, 0),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="domain-ge") == 60
    assert M.det(method="bareiss") == 60
    assert M.det(method="berkowitz") == 60
    assert M.det(method="lu") == 60

    M = Matrix(( ( 1,  0,  0,  0),
                 ( 5,  0,  0,  0),
                 ( 9, 10, 11, 0),
                 (13, 14, 15, 16) ))

    assert M.det(method="domain-ge") == 0
    assert M.det(method="bareiss") == 0
    assert M.det(method="berkowitz") == 0
    assert M.det(method="lu") == 0

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (0, 0, 0, 0, 3) ))

    assert M.det(method="domain-ge") == 243
    assert M.det(method="bareiss") == 243
    assert M.det(method="berkowitz") == 243
    assert M.det(method="lu") == 243

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="domain-ge") == -55
    assert M.det(method="bareiss") == -55
    assert M.det(method="berkowitz") == -55
    assert M.det(method="lu") == -55

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="domain-ge") == 11664
    assert M.det(method="bareiss") == 11664
    assert M.det(method="berkowitz") == 11664
    assert M.det(method="lu") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="domain-ge") == 123
    assert M.det(method="bareiss") == 123
    assert M.det(method="berkowitz") == 123
    assert M.det(method="lu") == 123

    M = Matrix(( (x, y, z),
                 (1, 0, 0),
                 (y, z, x) ))

    assert M.det(method="domain-ge") == z**2 - x*y
    assert M.det(method="bareiss") == z**2 - x*y
    assert M.det(method="berkowitz") == z**2 - x*y
    assert M.det(method="lu") == z**2 - x*y

    # issue 13835
    a = symbols('a')
    M = lambda n: Matrix([[i + a*j for i in range(n)]
                          for j in range(n)])
    assert M(5).det() == 0
    assert M(6).det() == 0
    assert M(7).det() == 0

def test_issue_14517():
    M = Matrix([
        [   0, 10*I,    10*I,       0],
        [10*I,    0,       0,    10*I],
        [10*I,    0, 5 + 2*I,    10*I],
        [   0, 10*I,    10*I, 5 + 2*I]])
    ev = M.eigenvals()
    # test one random eigenvalue, the computation is a little slow
    test_ev = random.choice(list(ev.keys()))
    assert (M - test_ev*eye(4)).det() == 0

def test_legacy_det():
    # Minimal support for legacy keys for 'method' in det()
    # Partially copied from test_determinant()

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="bareis") == -289
    assert M.det(method="det_lu") == -289
    assert M.det(method="det_LU") == -289

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="bareis") == 275
    assert M.det(method="det_lu") == 275
    assert M.det(method="Bareis") == 275

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="bareis") == -55
    assert M.det(method="det_lu") == -55
    assert M.det(method="BAREISS") == -55

    M = Matrix(( ( 3,  0,  0, 0),
                 (-2,  1,  0, 0),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="bareiss") == 60
    assert M.det(method="berkowitz") == 60
    assert M.det(method="lu") == 60

    M = Matrix(( ( 1,  0,  0,  0),
                 ( 5,  0,  0,  0),
                 ( 9, 10, 11, 0),
                 (13, 14, 15, 16) ))

    assert M.det(method="bareiss") == 0
    assert M.det(method="berkowitz") == 0
    assert M.det(method="lu") == 0

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (0, 0, 0, 0, 3) ))

    assert M.det(method="bareiss") == 243
    assert M.det(method="berkowitz") == 243
    assert M.det(method="lu") == 243

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="bareis") == 11664
    assert M.det(method="det_lu") == 11664
    assert M.det(method="BERKOWITZ") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="bareis") == 123
    assert M.det(method="det_lu") == 123
    assert M.det(method="LU") == 123

def eye_Determinant(n):
    return Matrix(n, n, lambda i, j: int(i == j))

def zeros_Determinant(n):
    return Matrix(n, n, lambda i, j: 0)

def test_det():
    a = Matrix(2, 3, [1, 2, 3, 4, 5, 6])
    raises(NonSquareMatrixError, lambda: a.det())

    z = zeros_Determinant(2)
    ey = eye_Determinant(2)
    assert z.det() == 0
    assert ey.det() == 1

    x = Symbol('x')
    a = Matrix(0, 0, [])
    b = Matrix(1, 1, [5])
    c = Matrix(2, 2, [1, 2, 3, 4])
    d = Matrix(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 8])
    e = Matrix(4, 4,
        [x, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])
    from sympy.abc import i, j, k, l, m, n
    f = Matrix(3, 3, [i, l, m, 0, j, n, 0, 0, k])
    g = Matrix(3, 3, [i, 0, 0, l, j, 0, m, n, k])
    h = Matrix(3, 3, [x**3, 0, 0, i, x**-1, 0, j, k, x**-2])
    # the method keyword for `det` doesn't kick in until 4x4 matrices,
    # so there is no need to test all methods on smaller ones

    assert a.det() == 1
    assert b.det() == 5
    assert c.det() == -2
    assert d.det() == 3
    assert e.det() == 4*x - 24
    assert e.det(method="domain-ge") == 4*x - 24
    assert e.det(method='bareiss') == 4*x - 24
    assert e.det(method='berkowitz') == 4*x - 24
    assert f.det() == i*j*k
    assert g.det() == i*j*k
    assert h.det() == 1
    raises(ValueError, lambda: e.det(iszerofunc="test"))

def test_permanent():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert M.per() == 450
    for i in range(1, 12):
        assert ones(i, i).per() == ones(i, i).T.per() == factorial(i)
        assert (ones(i, i)-eye(i)).per() == (ones(i, i)-eye(i)).T.per() == subfactorial(i)

    a1, a2, a3, a4, a5 = symbols('a_1 a_2 a_3 a_4 a_5')
    M = Matrix([a1, a2, a3, a4, a5])
    assert M.per() == M.T.per() == a1 + a2 + a3 + a4 + a5

def test_adjugate():
    x = Symbol('x')
    e = Matrix(4, 4,
        [x, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])

    adj = Matrix([
        [   4,         -8,         4,         0],
        [  76, -14*x - 68,  14*x - 8, -4*x + 24],
        [-122, 17*x + 142, -21*x + 4,  8*x - 48],
        [  48,  -4*x - 72,       8*x, -4*x + 24]])
    assert e.adjugate() == adj
    assert e.adjugate(method='bareiss') == adj
    assert e.adjugate(method='berkowitz') == adj

    a = Matrix(2, 3, [1, 2, 3, 4, 5, 6])
    raises(NonSquareMatrixError, lambda: a.adjugate())

def test_util():
    R = Rational

    v1 = Matrix(1, 3, [1, 2, 3])
    v2 = Matrix(1, 3, [3, 4, 5])
    assert v1.norm() == sqrt(14)
    assert v1.project(v2) == Matrix(1, 3, [R(39)/25, R(52)/25, R(13)/5])
    assert Matrix.zeros(1, 2) == Matrix(1, 2, [0, 0])
    assert ones(1, 2) == Matrix(1, 2, [1, 1])
    assert v1.copy() == v1
    # cofactor
    assert eye(3) == eye(3).cofactor_matrix()
    test = Matrix([[1, 3, 2], [2, 6, 3], [2, 3, 6]])
    assert test.cofactor_matrix() == \
        Matrix([[27, -6, -6], [-12, 2, 3], [-3, 1, 0]])
    test = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert test.cofactor_matrix() == \
        Matrix([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])

def test_cofactor_and_minors():
    x = Symbol('x')
    e = Matrix(4, 4,
        [x, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])

    m = Matrix([
        [ x,  1,  3],
        [ 2,  9, 11],
        [12, 13, 14]])
    cm = Matrix([
        [ 4,         76,       -122,        48],
        [-8, -14*x - 68, 17*x + 142, -4*x - 72],
        [ 4,   14*x - 8,  -21*x + 4,       8*x],
        [ 0,  -4*x + 24,   8*x - 48, -4*x + 24]])
    sub = Matrix([
            [x, 1,  2],
            [4, 5,  6],
            [2, 9, 10]])

    assert e.minor_submatrix(1, 2) == m
    assert e.minor_submatrix(-1, -1) == sub
    assert e.minor(1, 2) == -17*x - 142
    assert e.cofactor(1, 2) == 17*x + 142
    assert e.cofactor_matrix() == cm
    assert e.cofactor_matrix(method="bareiss") == cm
    assert e.cofactor_matrix(method="berkowitz") == cm

    raises(ValueError, lambda: e.cofactor(4, 5))
    raises(ValueError, lambda: e.minor(4, 5))
    raises(ValueError, lambda: e.minor_submatrix(4, 5))

    a = Matrix(2, 3, [1, 2, 3, 4, 5, 6])
    assert a.minor_submatrix(0, 0) == Matrix([[5, 6]])

    raises(ValueError, lambda:
        Matrix(0, 0, []).minor_submatrix(0, 0))
    raises(NonSquareMatrixError, lambda: a.cofactor(0, 0))
    raises(NonSquareMatrixError, lambda: a.minor(0, 0))
    raises(NonSquareMatrixError, lambda: a.cofactor_matrix())

def test_charpoly():
    x, y = Symbol('x'), Symbol('y')
    z, t = Symbol('z'), Symbol('t')

    from sympy.abc import a,b,c

    m = Matrix(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9])

    assert eye_Determinant(3).charpoly(x) == Poly((x - 1)**3, x)
    assert eye_Determinant(3).charpoly(y) == Poly((y - 1)**3, y)
    assert m.charpoly() == Poly(x**3 - 15*x**2 - 18*x, x)
    raises(NonSquareMatrixError, lambda: Matrix([[1], [2]]).charpoly())
    n = Matrix(4, 4, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert n.charpoly() == Poly(x**4, x)

    n = Matrix(4, 4, [45, 0, 0, 0, 0, 23, 0, 0, 0, 0, 87, 0, 0, 0, 0, 12])
    assert n.charpoly() == Poly(x**4 - 167*x**3 + 8811*x**2 - 173457*x + 1080540, x)

    n = Matrix(3, 3, [x, 0, 0, a, y, 0, b, c, z])
    assert n.charpoly() == Poly(t**3 - (x+y+z)*t**2 + t*(x*y+y*z+x*z) - x*y*z, t)
