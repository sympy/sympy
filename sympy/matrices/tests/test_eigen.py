from sympy import (
    Rational, Poly, Symbol, N, I, Abs, sqrt, exp, Float, sin,
    cos, symbols)
from sympy.matrices import eye, Matrix
from sympy.matrices.common import _MinimalMatrix, _CastableMatrix, MatrixEigen
from sympy.abc import x, y
from sympy.core.singleton import S
from sympy.testing.pytest import raises, XFAIL
from sympy.matrices.matrices import NonSquareMatrixError, MatrixError
from sympy.simplify.simplify import simplify
class EigenOnlyMatrix(_MinimalMatrix, _CastableMatrix, MatrixEigen):
    pass

def test_eigen():
    R = Rational

    assert eye(3).charpoly(x) == Poly((x - 1)**3, x)
    assert eye(3).charpoly(y) == Poly((y - 1)**3, y)

    M = Matrix([[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]])

    assert M.eigenvals(multiple=False) == {S.One: 3}
    assert M.eigenvals(multiple=True) == [1, 1, 1]

    assert M.eigenvects() == (
        [(1, 3, [Matrix([1, 0, 0]),
                 Matrix([0, 1, 0]),
                 Matrix([0, 0, 1])])])

    assert M.left_eigenvects() == (
        [(1, 3, [Matrix([[1, 0, 0]]),
                 Matrix([[0, 1, 0]]),
                 Matrix([[0, 0, 1]])])])

    M = Matrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])

    assert M.eigenvals() == {2*S.One: 1, -S.One: 1, S.Zero: 1}

    assert M.eigenvects() == (
        [
            (-1, 1, [Matrix([-1, 1, 0])]),
            ( 0, 1, [Matrix([0, -1, 1])]),
            ( 2, 1, [Matrix([R(2, 3), R(1, 3), 1])])
        ])

    assert M.left_eigenvects() == (
        [
            (-1, 1, [Matrix([[-2, 1, 1]])]),
            (0, 1, [Matrix([[-1, -1, 1]])]),
            (2, 1, [Matrix([[1, 1, 1]])])
        ])

    a = Symbol('a')
    M = Matrix([[a, 0],
                [0, 1]])

    assert M.eigenvals() == {a: 1, S.One: 1}

    M = Matrix([[1, -1],
                [1,  3]])
    assert M.eigenvects() == ([(2, 2, [Matrix(2, 1, [-1, 1])])])
    assert M.left_eigenvects() == ([(2, 2, [Matrix([[1, 1]])])])

    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    a = R(15, 2)
    b = 3*33**R(1, 2)
    c = R(13, 2)
    d = (R(33, 8) + 3*b/8)
    e = (R(33, 8) - 3*b/8)

    def NS(e, n):
        return str(N(e, n))
    r = [
        (a - b/2, 1, [Matrix([(12 + 24/(c - b/2))/((c - b/2)*e) + 3/(c - b/2),
                              (6 + 12/(c - b/2))/e, 1])]),
        (      0, 1, [Matrix([1, -2, 1])]),
        (a + b/2, 1, [Matrix([(12 + 24/(c + b/2))/((c + b/2)*d) + 3/(c + b/2),
                              (6 + 12/(c + b/2))/d, 1])]),
    ]
    r1 = [(NS(r[i][0], 2), NS(r[i][1], 2),
        [NS(j, 2) for j in r[i][2][0]]) for i in range(len(r))]
    r = M.eigenvects()
    r2 = [(NS(r[i][0], 2), NS(r[i][1], 2),
        [NS(j, 2) for j in r[i][2][0]]) for i in range(len(r))]
    assert sorted(r1) == sorted(r2)

    eps = Symbol('eps', real=True)

    M = Matrix([[abs(eps), I*eps    ],
                [-I*eps,   abs(eps) ]])

    assert M.eigenvects() == (
        [
            ( 0, 1, [Matrix([[-I*eps/abs(eps)], [1]])]),
            ( 2*abs(eps), 1, [ Matrix([[I*eps/abs(eps)], [1]]) ] ),
        ])

    assert M.left_eigenvects() == (
        [
            (0, 1, [Matrix([[I*eps/Abs(eps), 1]])]),
            (2*Abs(eps), 1, [Matrix([[-I*eps/Abs(eps), 1]])])
        ])

    M = Matrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
    M._eigenvects = M.eigenvects(simplify=False)
    assert max(i.q for i in M._eigenvects[0][2][0]) > 1
    M._eigenvects = M.eigenvects(simplify=True)
    assert max(i.q for i in M._eigenvects[0][2][0]) == 1
    M = Matrix([[Rational(1, 4), 1], [1, 1]])
    assert M.eigenvects(simplify=True) == [
        (Rational(5, 8) - sqrt(73)/8, 1, [Matrix([[-sqrt(73)/8 - Rational(3, 8)], [1]])]),
        (Rational(5, 8) + sqrt(73)/8, 1, [Matrix([[Rational(-3, 8) + sqrt(73)/8], [1]])])]
    assert M.eigenvects(simplify=False) == [
        (Rational(5, 8) - sqrt(73)/8, 1, [Matrix([[-1/(-Rational(3, 8) + sqrt(73)/8)], [1]])]),
        (Rational(5, 8) + sqrt(73)/8, 1, [Matrix([[8/(3 + sqrt(73))], [1]])])]

    m = Matrix([[1, .6, .6], [.6, .9, .9], [.9, .6, .6]])
    evals = { Rational(5, 4) - sqrt(385)/20: 1, sqrt(385)/20 + Rational(5, 4): 1, S.Zero: 1}
    assert m.eigenvals() == evals
    nevals = list(sorted(m.eigenvals(rational=False).keys()))
    sevals = list(sorted(evals.keys()))
    assert all(abs(nevals[i] - sevals[i]) < 1e-9 for i in range(len(nevals)))

    # issue 10719
    assert Matrix([]).eigenvals() == {}
    assert Matrix([]).eigenvects() == []

    # issue 15119
    raises(NonSquareMatrixError, lambda : Matrix([[1, 2], [0, 4], [0, 0]]).eigenvals())
    raises(NonSquareMatrixError, lambda : Matrix([[1, 0], [3, 4], [5, 6]]).eigenvals())
    raises(NonSquareMatrixError, lambda : Matrix([[1, 2, 3], [0, 5, 6]]).eigenvals())
    raises(NonSquareMatrixError, lambda : Matrix([[1, 0, 0], [4, 5, 0]]).eigenvals())
    raises(NonSquareMatrixError, lambda : Matrix([[1, 2, 3], [0, 5, 6]]).eigenvals(error_when_incomplete = False))
    raises(NonSquareMatrixError, lambda : Matrix([[1, 0, 0], [4, 5, 0]]).eigenvals(error_when_incomplete = False))

    # issue 15125
    from sympy.core.function import count_ops
    q = Symbol("q", positive = True)
    m = Matrix([[-2, exp(-q), 1], [exp(q), -2, 1], [1, 1, -2]])
    assert count_ops(m.eigenvals(simplify=False)) > count_ops(m.eigenvals(simplify=True))
    assert count_ops(m.eigenvals(simplify=lambda x: x)) > count_ops(m.eigenvals(simplify=True))

    assert isinstance(m.eigenvals(simplify=True, multiple=False), dict)
    assert isinstance(m.eigenvals(simplify=True, multiple=True), list)
    assert isinstance(m.eigenvals(simplify=lambda x: x, multiple=False), dict)
    assert isinstance(m.eigenvals(simplify=lambda x: x, multiple=True), list)

@XFAIL
def test_eigen_vects():
    m = Matrix(2, 2, [1, 0, 0, I])
    raises(NotImplementedError, lambda: m.is_diagonalizable(True))
    # !!! bug because of eigenvects() or roots(x**2 + (-1 - I)*x + I, x)
    # see issue 5292
    assert not m.is_diagonalizable(True)
    raises(MatrixError, lambda: m.diagonalize(True))
    (P, D) = m.diagonalize(True)

def test_issue_8240():
    # Eigenvalues of large triangular matrices
    n = 200

    diagonal_variables = [Symbol('x%s' % i) for i in range(n)]
    M = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        M[i][i] = diagonal_variables[i]
    M = Matrix(M)

    eigenvals = M.eigenvals()
    assert len(eigenvals) == n
    for i in range(n):
        assert eigenvals[diagonal_variables[i]] == 1

    eigenvals = M.eigenvals(multiple=True)
    assert set(eigenvals) == set(diagonal_variables)

    # with multiplicity
    M = Matrix([[x, 0, 0], [1, y, 0], [2, 3, x]])
    eigenvals = M.eigenvals()
    assert eigenvals == {x: 2, y: 1}

    eigenvals = M.eigenvals(multiple=True)
    assert len(eigenvals) == 3
    assert eigenvals.count(x) == 2
    assert eigenvals.count(y) == 1

# EigenOnlyMatrix tests
def test_eigenvals():
    M = EigenOnlyMatrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])
    assert M.eigenvals() == {2*S.One: 1, -S.One: 1, S.Zero: 1}

    # if we cannot factor the char poly, we raise an error
    m = Matrix([
        [3,  0,  0, 0, -3],
        [0, -3, -3, 0,  3],
        [0,  3,  0, 3,  0],
        [0,  0,  3, 0,  3],
        [3,  0,  0, 3,  0]])
    raises(MatrixError, lambda: m.eigenvals())


def test_eigenvects():
    M = EigenOnlyMatrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])
    vecs = M.eigenvects()
    for val, mult, vec_list in vecs:
        assert len(vec_list) == 1
        assert M*vec_list[0] == val*vec_list[0]


def test_left_eigenvects():
    M = EigenOnlyMatrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])
    vecs = M.left_eigenvects()
    for val, mult, vec_list in vecs:
        assert len(vec_list) == 1
        assert vec_list[0]*M == val*vec_list[0]


def test_diagonalize():
    m = EigenOnlyMatrix(2, 2, [0, -1, 1, 0])
    raises(MatrixError, lambda: m.diagonalize(reals_only=True))
    P, D = m.diagonalize()
    assert D.is_diagonal()
    assert D == Matrix([
                 [-I, 0],
                 [ 0, I]])

    # make sure we use floats out if floats are passed in
    m = EigenOnlyMatrix(2, 2, [0, .5, .5, 0])
    P, D = m.diagonalize()
    assert all(isinstance(e, Float) for e in D.values())
    assert all(isinstance(e, Float) for e in P.values())

    _, D2 = m.diagonalize(reals_only=True)
    assert D == D2


def test_is_diagonalizable():
    a, b, c = symbols('a b c')
    m = EigenOnlyMatrix(2, 2, [a, c, c, b])
    assert m.is_symmetric()
    assert m.is_diagonalizable()
    assert not EigenOnlyMatrix(2, 2, [1, 1, 0, 1]).is_diagonalizable()

    m = EigenOnlyMatrix(2, 2, [0, -1, 1, 0])
    assert m.is_diagonalizable()
    assert not m.is_diagonalizable(reals_only=True)


def test_jordan_form():
    m = Matrix(3, 2, [-3, 1, -3, 20, 3, 10])
    raises(NonSquareMatrixError, lambda: m.jordan_form())

    # the next two tests test the cases where the old
    # algorithm failed due to the fact that the block structure can
    # *NOT* be determined  from algebraic and geometric multiplicity alone
    # This can be seen most easily when one lets compute the J.c.f. of a matrix that
    # is in J.c.f already.
    m = EigenOnlyMatrix(4, 4, [2, 1, 0, 0,
                    0, 2, 1, 0,
                    0, 0, 2, 0,
                    0, 0, 0, 2
    ])
    P, J = m.jordan_form()
    assert m == J

    m = EigenOnlyMatrix(4, 4, [2, 1, 0, 0,
                    0, 2, 0, 0,
                    0, 0, 2, 1,
                    0, 0, 0, 2
    ])
    P, J = m.jordan_form()
    assert m == J

    A = Matrix([[ 2,  4,  1,  0],
                [-4,  2,  0,  1],
                [ 0,  0,  2,  4],
                [ 0,  0, -4,  2]])
    P, J = A.jordan_form()
    assert simplify(P*J*P.inv()) == A

    assert EigenOnlyMatrix(1, 1, [1]).jordan_form() == (
        Matrix([1]), Matrix([1]))
    assert EigenOnlyMatrix(1, 1, [1]).jordan_form(
        calc_transform=False) == Matrix([1])

    # make sure if we cannot factor the characteristic polynomial, we raise an error
    m = Matrix([[3, 0, 0, 0, -3], [0, -3, -3, 0, 3], [0, 3, 0, 3, 0], [0, 0, 3, 0, 3], [3, 0, 0, 3, 0]])
    raises(MatrixError, lambda: m.jordan_form())

    # make sure that if the input has floats, the output does too
    m = Matrix([
        [                0.6875, 0.125 + 0.1875*sqrt(3)],
        [0.125 + 0.1875*sqrt(3),                 0.3125]])
    P, J = m.jordan_form()
    assert all(isinstance(x, Float) or x == 0 for x in P)
    assert all(isinstance(x, Float) or x == 0 for x in J)


def test_singular_values():
    x = Symbol('x', real=True)

    A = EigenOnlyMatrix([[0, 1*I], [2, 0]])
    # if singular values can be sorted, they should be in decreasing order
    assert A.singular_values() == [2, 1]

    A = eye(3)
    A[1, 1] = x
    A[2, 2] = 5
    vals = A.singular_values()
    # since Abs(x) cannot be sorted, test set equality
    assert set(vals) == set([5, 1, Abs(x)])

    A = EigenOnlyMatrix([[sin(x), cos(x)], [-cos(x), sin(x)]])
    vals = [sv.trigsimp() for sv in A.singular_values()]
    assert vals == [S.One, S.One]

    A = EigenOnlyMatrix([
        [2, 4],
        [1, 3],
        [0, 0],
        [0, 0]
        ])
    assert A.singular_values() == \
        [sqrt(sqrt(221) + 15), sqrt(15 - sqrt(221))]
    assert A.T.singular_values() == \
        [sqrt(sqrt(221) + 15), sqrt(15 - sqrt(221)), 0, 0]

def test___eq__():
    assert (EigenOnlyMatrix(
        [[0, 1, 1],
        [1, 0, 0],
        [1, 1, 1]]) == {}) is False
