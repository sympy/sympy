from sympy import MatrixSymbol, Q, ask, Identity, ZeroMatrix

X = MatrixSymbol('X', 2, 2)
Y = MatrixSymbol('Y', 2, 3)
Z = MatrixSymbol('Z', 2, 2)

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

def test_symmetric():
    assert ask(Q.symmetric(X), Q.symmetric(X))
    assert ask(Q.symmetric(X*Z), Q.symmetric(X)) is None
    assert ask(Q.symmetric(X*Z), Q.symmetric(X) & Q.symmetric(Z)) is True
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
            Q.positive_definite(X) & Q.orthogonal(Y)) is True
    assert ask(Q.positive_definite(Identity(3))) is True
    assert ask(Q.positive_definite(ZeroMatrix(3, 3))) is False
