from sympy import MatrixSymbol, Q, ask, Identity, ZeroMatrix, Trace, MatrixSlice
from sympy.utilities.pytest import XFAIL
from sympy.assumptions import assuming

X = MatrixSymbol('X', 2, 2)
Y = MatrixSymbol('Y', 2, 3)
Z = MatrixSymbol('Z', 2, 2)

def test_square():
    assert ask(Q.square(X))
    assert not ask(Q.square(Y))
    assert ask(Q.square(Y*Y.T))

def test_invertible():
    assert ask(Q.invertible(X), Q.invertible(X))
    assert ask(Q.invertible(Y)) is False
    assert ask(Q.invertible(X*Y), Q.invertible(X)) is False
    assert ask(Q.invertible(X*Z), Q.invertible(X)) is None
    assert ask(Q.invertible(X*Z), Q.invertible(X) & Q.invertible(Z)) is True
    assert ask(Q.invertible(X.T)) is None
    assert ask(Q.invertible(X.T), Q.invertible(X)) is True
    assert ask(Q.invertible(X.I)) is True
    assert ask(Q.invertible(Identity(3))) is True
    assert ask(Q.invertible(ZeroMatrix(3, 3))) is False
    assert ask(Q.invertible(X), Q.fullrank(X) & Q.square(X))

def test_singular():
    assert ask(Q.singular(X)) is None
    assert ask(Q.singular(X), Q.invertible(X)) is False
    assert ask(Q.singular(X), ~Q.invertible(X)) is True

@XFAIL
def test_invertible_fullrank():
    assert ask(Q.invertible(X), Q.fullrank(X))


def test_symmetric():
    assert ask(Q.symmetric(X), Q.symmetric(X))
    assert ask(Q.symmetric(X*Z), Q.symmetric(X)) is None
    assert ask(Q.symmetric(X*Z), Q.symmetric(X) & Q.symmetric(Z)) is True
    assert ask(Q.symmetric(X + Z), Q.symmetric(X) & Q.symmetric(Z)) is True
    assert ask(Q.symmetric(Y)) is False
    assert ask(Q.symmetric(Y*Y.T)) is True
    assert ask(Q.symmetric(Y.T*X*Y)) is None
    assert ask(Q.symmetric(Y.T*X*Y), Q.symmetric(X)) is True
    assert ask(Q.symmetric(X*X*X*X*X*X*X*X*X*X), Q.symmetric(X)) is True


def test_orthogonal():
    assert ask(Q.orthogonal(X), Q.orthogonal(X))
    assert ask(Q.orthogonal(X.T), Q.orthogonal(X)) is True
    assert ask(Q.orthogonal(X.I), Q.orthogonal(X)) is True
    assert ask(Q.orthogonal(Y)) is False
    assert ask(Q.orthogonal(X)) is None
    assert ask(Q.orthogonal(X*Z*X), Q.orthogonal(X) & Q.orthogonal(Z)) is True
    assert ask(Q.orthogonal(Identity(3))) is True
    assert ask(Q.orthogonal(ZeroMatrix(3, 3))) is False
    assert ask(Q.invertible(X), Q.orthogonal(X))
    assert not ask(Q.orthogonal(X + Z), Q.orthogonal(X) & Q.orthogonal(Z))

def test_fullrank():
    assert ask(Q.fullrank(X), Q.fullrank(X))
    assert ask(Q.fullrank(X.T), Q.fullrank(X)) is True
    assert ask(Q.fullrank(X)) is None
    assert ask(Q.fullrank(Y)) is None
    assert ask(Q.fullrank(X*Z), Q.fullrank(X) & Q.fullrank(Z)) is True
    assert ask(Q.fullrank(Identity(3))) is True
    assert ask(Q.fullrank(ZeroMatrix(3, 3))) is False
    assert ask(Q.invertible(X), ~Q.fullrank(X)) == False


def test_positive_definite():
    assert ask(Q.positive_definite(X), Q.positive_definite(X))
    assert ask(Q.positive_definite(X.T), Q.positive_definite(X)) is True
    assert ask(Q.positive_definite(X.I), Q.positive_definite(X)) is True
    assert ask(Q.positive_definite(Y)) is False
    assert ask(Q.positive_definite(X)) is None
    assert ask(Q.positive_definite(X*Z*X),
            Q.positive_definite(X) & Q.positive_definite(Z)) is True
    assert ask(Q.positive_definite(X), Q.orthogonal(X))
    assert ask(Q.positive_definite(Y.T*X*Y),
            Q.positive_definite(X) & Q.fullrank(Y)) is True
    assert not ask(Q.positive_definite(Y.T*X*Y), Q.positive_definite(X))
    assert ask(Q.positive_definite(Identity(3))) is True
    assert ask(Q.positive_definite(ZeroMatrix(3, 3))) is False
    assert ask(Q.positive_definite(X + Z), Q.positive_definite(X) &
            Q.positive_definite(Z)) is True
    assert not ask(Q.positive_definite(-X), Q.positive_definite(X))


def test_triangular():
    assert ask(Q.upper_triangular(X + Z.T + Identity(2)), Q.upper_triangular(X) &
            Q.lower_triangular(Z)) is True
    assert ask(Q.upper_triangular(X*Z.T), Q.upper_triangular(X) &
            Q.lower_triangular(Z)) is True
    assert ask(Q.lower_triangular(Identity(3))) is True
    assert ask(Q.lower_triangular(ZeroMatrix(3, 3))) is True
    assert ask(Q.triangular(X), Q.unit_triangular(X))


def test_diagonal():
    assert ask(Q.diagonal(X + Z.T + Identity(2)), Q.diagonal(X) &
               Q.diagonal(Z)) is True
    assert ask(Q.diagonal(ZeroMatrix(3, 3)))
    assert ask(Q.lower_triangular(X) & Q.upper_triangular(X), Q.diagonal(X))
    assert ask(Q.diagonal(X), Q.lower_triangular(X) & Q.upper_triangular(X))
    assert ask(Q.symmetric(X), Q.diagonal(X))
    assert ask(Q.triangular(X), Q.diagonal(X))


def test_non_atoms():
    assert ask(Q.real(Trace(X)), Q.positive(Trace(X)))

@XFAIL
def test_non_trivial_implies():
    X = MatrixSymbol('X', 3, 3)
    Y = MatrixSymbol('Y', 3, 3)
    assert ask(Q.lower_triangular(X+Y), Q.lower_triangular(X) &
               Q.lower_triangular(Y))
    assert ask(Q.triangular(X), Q.lower_triangular(X))
    assert ask(Q.triangular(X+Y), Q.lower_triangular(X) &
               Q.lower_triangular(Y))

def test_MatrixSlice():
    X = MatrixSymbol('X', 4, 4)
    B = MatrixSlice(X, (1, 3), (1, 3))
    C = MatrixSlice(X, (0, 3), (1, 3))
    assert ask(Q.symmetric(B), Q.symmetric(X))
    assert ask(Q.invertible(B), Q.invertible(X))
    assert ask(Q.diagonal(B), Q.diagonal(X))
    assert ask(Q.orthogonal(B), Q.orthogonal(X))
    assert ask(Q.upper_triangular(B), Q.upper_triangular(X))

    assert not ask(Q.symmetric(C), Q.symmetric(X))
    assert not ask(Q.invertible(C), Q.invertible(X))
    assert not ask(Q.diagonal(C), Q.diagonal(X))
    assert not ask(Q.orthogonal(C), Q.orthogonal(X))
    assert not ask(Q.upper_triangular(C), Q.upper_triangular(X))
