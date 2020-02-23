from sympy.matrices.common import _MinimalMatrix, _CastableMatrix
from sympy.matrices.matrices import MatrixSubspaces
from sympy.matrices import Matrix

class SubspaceOnlyMatrix(_MinimalMatrix, _CastableMatrix, MatrixSubspaces):
    pass


def test_columnspace():
    m = SubspaceOnlyMatrix([[ 1,  2,  0,  2,  5],
                            [-2, -5,  1, -1, -8],
                            [ 0, -3,  3,  4,  1],
                            [ 3,  6,  0, -7,  2]])

    basis = m.columnspace()
    assert basis[0] == Matrix([1, -2, 0, 3])
    assert basis[1] == Matrix([2, -5, -3, 6])
    assert basis[2] == Matrix([2, -1, 4, -7])

    assert len(basis) == 3
    assert Matrix.hstack(m, *basis).columnspace() == basis


def test_rowspace():
    m = SubspaceOnlyMatrix([[ 1,  2,  0,  2,  5],
                            [-2, -5,  1, -1, -8],
                            [ 0, -3,  3,  4,  1],
                            [ 3,  6,  0, -7,  2]])

    basis = m.rowspace()
    assert basis[0] == Matrix([[1, 2, 0, 2, 5]])
    assert basis[1] == Matrix([[0, -1, 1, 3, 2]])
    assert basis[2] == Matrix([[0, 0, 0, 5, 5]])

    assert len(basis) == 3


def test_nullspace():
    m = SubspaceOnlyMatrix([[ 1,  2,  0,  2,  5],
                            [-2, -5,  1, -1, -8],
                            [ 0, -3,  3,  4,  1],
                            [ 3,  6,  0, -7,  2]])

    basis = m.nullspace()
    assert basis[0] == Matrix([-2, 1, 1, 0, 0])
    assert basis[1] == Matrix([-1, -1, 0, -1, 1])
    # make sure the null space is really gets zeroed
    assert all(e.is_zero for e in m*basis[0])
    assert all(e.is_zero for e in m*basis[1])


def test_orthogonalize():
    m = Matrix([[1, 2], [3, 4]])
    assert m.orthogonalize(Matrix([[2], [1]])) == [Matrix([[2], [1]])]
    assert m.orthogonalize(Matrix([[2], [1]]), normalize=True) == \
        [Matrix([[2*sqrt(5)/5], [sqrt(5)/5]])]
    assert m.orthogonalize(Matrix([[1], [2]]), Matrix([[-1], [4]])) == \
        [Matrix([[1], [2]]), Matrix([[Rational(-12, 5)], [Rational(6, 5)]])]
    assert m.orthogonalize(Matrix([[0], [0]]), Matrix([[-1], [4]])) == \
        [Matrix([[-1], [4]])]
    assert m.orthogonalize(Matrix([[0], [0]])) == []

    n = Matrix([[9, 1, 9], [3, 6, 10], [8, 5, 2]])
    vecs = [Matrix([[-5], [1]]), Matrix([[-5], [2]]), Matrix([[-5], [-2]])]
    assert n.orthogonalize(*vecs) == \
        [Matrix([[-5], [1]]), Matrix([[Rational(5, 26)], [Rational(25, 26)]])]

    vecs = [Matrix([0, 0, 0]), Matrix([1, 2, 3]), Matrix([1, 4, 5])]
    raises(ValueError, lambda: Matrix.orthogonalize(*vecs, rankcheck=True))

    vecs = [Matrix([1, 2, 3]), Matrix([4, 5, 6]), Matrix([7, 8, 9])]
    raises(ValueError, lambda: Matrix.orthogonalize(*vecs, rankcheck=True))
