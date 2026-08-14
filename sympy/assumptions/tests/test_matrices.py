from __future__ import annotations
from sympy.assumptions.ask import (Q, _ask_recursive)
from sympy.core.symbol import Symbol
from sympy.matrices.expressions.diagonal import (DiagMatrix, DiagonalMatrix)
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions import (MatrixSymbol, Identity, ZeroMatrix,
        OneMatrix, Trace, MatrixSlice, Determinant, BlockMatrix, BlockDiagMatrix)
from sympy.matrices.expressions.factorizations import LofLU
from sympy.testing.pytest import XFAIL

X = MatrixSymbol('X', 2, 2)
Y = MatrixSymbol('Y', 2, 3)
Z = MatrixSymbol('Z', 2, 2)
A1x1 = MatrixSymbol('A1x1', 1, 1)
B1x1 = MatrixSymbol('B1x1', 1, 1)
C0x0 = MatrixSymbol('C0x0', 0, 0)
V1 = MatrixSymbol('V1', 2, 1)
V2 = MatrixSymbol('V2', 2, 1)

def test_square():
    assert _ask_recursive(Q.square(X)) is True
    assert _ask_recursive(Q.square(Y)) is False
    assert _ask_recursive(Q.square(Y*Y.T)) is True

def test_invertible():
    assert _ask_recursive(Q.invertible(X), Q.invertible(X)) is True
    assert _ask_recursive(Q.invertible(Y)) is False
    assert _ask_recursive(Q.invertible(X*Y), Q.invertible(X)) is False
    assert _ask_recursive(Q.invertible(X*Z), Q.invertible(X)) is None
    assert _ask_recursive(Q.invertible(X*Z), Q.invertible(X) & Q.invertible(Z)) is True
    assert _ask_recursive(Q.invertible(X.T)) is None
    assert _ask_recursive(Q.invertible(X.T), Q.invertible(X)) is True
    assert _ask_recursive(Q.invertible(X.I)) is True
    assert _ask_recursive(Q.invertible(Identity(3))) is True
    assert _ask_recursive(Q.invertible(ZeroMatrix(3, 3))) is False
    assert _ask_recursive(Q.invertible(OneMatrix(1, 1))) is True
    assert _ask_recursive(Q.invertible(OneMatrix(3, 3))) is False
    assert _ask_recursive(Q.invertible(X), Q.fullrank(X) & Q.square(X)) is None

def test_singular():
    assert _ask_recursive(Q.singular(X)) is None
    assert _ask_recursive(Q.singular(X), Q.invertible(X)) is False
    assert _ask_recursive(Q.singular(X), ~Q.invertible(X)) is None

@XFAIL
def test_invertible_fullrank():
    assert _ask_recursive(Q.invertible(X), Q.fullrank(X)) is True


def test_invertible_BlockMatrix():
    assert _ask_recursive(Q.invertible(BlockMatrix([Identity(3)]))) == True
    assert _ask_recursive(Q.invertible(BlockMatrix([ZeroMatrix(3, 3)]))) == False

    X = Matrix([[1, 2, 3], [3, 5, 4]])
    Y = Matrix([[4, 2, 7], [2, 3, 5]])
    # non-invertible A block
    assert _ask_recursive(Q.invertible(BlockMatrix([
        [Matrix.ones(3, 3), Y.T],
        [X, Matrix.eye(2)],
    ]))) == True
    # non-invertible B block
    assert _ask_recursive(Q.invertible(BlockMatrix([
        [Y.T, Matrix.ones(3, 3)],
        [Matrix.eye(2), X],
    ]))) == True
    # non-invertible C block
    assert _ask_recursive(Q.invertible(BlockMatrix([
        [X, Matrix.eye(2)],
        [Matrix.ones(3, 3), Y.T],
    ]))) == True
    # non-invertible D block
    assert _ask_recursive(Q.invertible(BlockMatrix([
        [Matrix.eye(2), X],
        [Y.T, Matrix.ones(3, 3)],
    ]))) == True


def test_invertible_BlockDiagMatrix():
    assert _ask_recursive(Q.invertible(BlockDiagMatrix(Identity(3), Identity(5)))) == True
    assert _ask_recursive(Q.invertible(BlockDiagMatrix(ZeroMatrix(3, 3), Identity(5)))) == False
    assert _ask_recursive(Q.invertible(BlockDiagMatrix(Identity(3), OneMatrix(5, 5)))) == False


def test_symmetric():
    assert _ask_recursive(Q.symmetric(X), Q.symmetric(X)) is True
    assert _ask_recursive(Q.symmetric(X*Z), Q.symmetric(X)) is None
    assert _ask_recursive(Q.symmetric(X*Z), Q.symmetric(X) & Q.symmetric(Z)) is True
    assert _ask_recursive(Q.symmetric(X + Z), Q.symmetric(X) & Q.symmetric(Z)) is True
    assert _ask_recursive(Q.symmetric(Y)) is False
    assert _ask_recursive(Q.symmetric(Y*Y.T)) is True
    assert _ask_recursive(Q.symmetric(Y.T*X*Y)) is None
    assert _ask_recursive(Q.symmetric(Y.T*X*Y), Q.symmetric(X)) is True
    assert _ask_recursive(Q.symmetric(X**10), Q.symmetric(X)) is True
    assert _ask_recursive(Q.symmetric(A1x1)) is True
    assert _ask_recursive(Q.symmetric(A1x1 + B1x1)) is True
    assert _ask_recursive(Q.symmetric(A1x1 * B1x1)) is True
    assert _ask_recursive(Q.symmetric(V1.T*V1)) is True
    assert _ask_recursive(Q.symmetric(V1.T*(V1 + V2))) is True
    assert _ask_recursive(Q.symmetric(V1.T*(V1 + V2) + A1x1)) is True
    assert _ask_recursive(Q.symmetric(MatrixSlice(Y, (0, 1), (1, 2)))) is True
    assert _ask_recursive(Q.symmetric(Identity(3))) is True
    assert _ask_recursive(Q.symmetric(ZeroMatrix(3, 3))) is True
    assert _ask_recursive(Q.symmetric(OneMatrix(3, 3))) is True

def _test_orthogonal_unitary(predicate):
    assert _ask_recursive(predicate(X), predicate(X)) is True
    assert _ask_recursive(predicate(X.T), predicate(X)) is True
    assert _ask_recursive(predicate(X.I), predicate(X)) is True
    assert _ask_recursive(predicate(X**2), predicate(X)) is True
    assert _ask_recursive(predicate(Y)) is False
    assert _ask_recursive(predicate(X)) is None
    assert _ask_recursive(predicate(X), ~Q.invertible(X)) is False
    assert _ask_recursive(predicate(X*Z*X), predicate(X) & predicate(Z)) is True
    assert _ask_recursive(predicate(Identity(3))) is True
    assert _ask_recursive(predicate(ZeroMatrix(3, 3))) is False
    assert _ask_recursive(Q.invertible(X), predicate(X)) is True
    assert _ask_recursive(predicate(X + Z), predicate(X) & predicate(Z)) is None

def test_orthogonal():
    _test_orthogonal_unitary(Q.orthogonal)

def test_unitary():
    _test_orthogonal_unitary(Q.unitary)
    assert _ask_recursive(Q.unitary(X), Q.orthogonal(X)) is True

def test_fullrank():
    assert _ask_recursive(Q.fullrank(X), Q.fullrank(X)) is True
    assert _ask_recursive(Q.fullrank(X**2), Q.fullrank(X)) is True
    assert _ask_recursive(Q.fullrank(X.T), Q.fullrank(X)) is True
    assert _ask_recursive(Q.fullrank(X)) is None
    assert _ask_recursive(Q.fullrank(Y)) is None
    assert _ask_recursive(Q.fullrank(X*Z), Q.fullrank(X) & Q.fullrank(Z)) is True
    assert _ask_recursive(Q.fullrank(Identity(3))) is True
    assert _ask_recursive(Q.fullrank(ZeroMatrix(3, 3))) is False
    assert _ask_recursive(Q.fullrank(OneMatrix(1, 1))) is True
    assert _ask_recursive(Q.fullrank(OneMatrix(3, 3))) is False
    assert _ask_recursive(Q.invertible(X), ~Q.fullrank(X)) == False


def test_positive_definite():
    assert _ask_recursive(Q.positive_definite(X), Q.positive_definite(X)) is True
    assert _ask_recursive(Q.positive_definite(X.T), Q.positive_definite(X)) is True
    assert _ask_recursive(Q.positive_definite(X.I), Q.positive_definite(X)) is True
    assert _ask_recursive(Q.positive_definite(Y)) is False
    assert _ask_recursive(Q.positive_definite(X)) is None
    assert _ask_recursive(Q.positive_definite(X**3), Q.positive_definite(X)) is True
    assert _ask_recursive(Q.positive_definite(X*Z*X),
            Q.positive_definite(X) & Q.positive_definite(Z)) is True
    assert _ask_recursive(Q.positive_definite(X), Q.orthogonal(X)) is True
    assert _ask_recursive(Q.positive_definite(Y.T*X*Y),
            Q.positive_definite(X) & Q.fullrank(Y)) is True
    assert _ask_recursive(Q.positive_definite(Y.T*X*Y), Q.positive_definite(X)) is None
    assert _ask_recursive(Q.positive_definite(Identity(3))) is True
    assert _ask_recursive(Q.positive_definite(ZeroMatrix(3, 3))) is False
    assert _ask_recursive(Q.positive_definite(OneMatrix(1, 1))) is True
    assert _ask_recursive(Q.positive_definite(OneMatrix(3, 3))) is False
    assert _ask_recursive(Q.positive_definite(X + Z), Q.positive_definite(X) &
            Q.positive_definite(Z)) is True
    assert _ask_recursive(Q.positive_definite(-X), Q.positive_definite(X)) is None
    assert _ask_recursive(Q.positive(X[1, 1]), Q.positive_definite(X)) is True

def test_triangular():
    assert _ask_recursive(Q.upper_triangular(X + Z.T + Identity(2)), Q.upper_triangular(X) &
            Q.lower_triangular(Z)) is True
    assert _ask_recursive(Q.upper_triangular(X*Z.T), Q.upper_triangular(X) &
            Q.lower_triangular(Z)) is True
    assert _ask_recursive(Q.lower_triangular(Identity(3))) is True
    assert _ask_recursive(Q.lower_triangular(ZeroMatrix(3, 3))) is True
    assert _ask_recursive(Q.upper_triangular(ZeroMatrix(3, 3))) is True
    assert _ask_recursive(Q.lower_triangular(OneMatrix(1, 1))) is True
    assert _ask_recursive(Q.upper_triangular(OneMatrix(1, 1))) is True
    assert _ask_recursive(Q.lower_triangular(OneMatrix(3, 3))) is False
    assert _ask_recursive(Q.upper_triangular(OneMatrix(3, 3))) is False
    assert _ask_recursive(Q.triangular(X), Q.unit_triangular(X)) is True
    assert _ask_recursive(Q.upper_triangular(X**3), Q.upper_triangular(X)) is True
    assert _ask_recursive(Q.lower_triangular(X**3), Q.lower_triangular(X)) is True


def test_diagonal():
    assert _ask_recursive(Q.diagonal(X + Z.T + Identity(2)), Q.diagonal(X) &
               Q.diagonal(Z)) is True
    assert _ask_recursive(Q.diagonal(ZeroMatrix(3, 3))) is True
    assert _ask_recursive(Q.diagonal(OneMatrix(1, 1))) is True
    assert _ask_recursive(Q.diagonal(OneMatrix(3, 3))) is False
    assert _ask_recursive(Q.lower_triangular(X) & Q.upper_triangular(X), Q.diagonal(X)) is True
    assert _ask_recursive(Q.diagonal(X), Q.lower_triangular(X) & Q.upper_triangular(X)) is None
    assert _ask_recursive(Q.symmetric(X), Q.diagonal(X)) is True
    assert _ask_recursive(Q.triangular(X), Q.diagonal(X)) is True
    assert _ask_recursive(Q.diagonal(C0x0)) is True
    assert _ask_recursive(Q.diagonal(A1x1)) is True
    assert _ask_recursive(Q.diagonal(A1x1 + B1x1)) is True
    assert _ask_recursive(Q.diagonal(A1x1*B1x1)) is True
    assert _ask_recursive(Q.diagonal(V1.T*V2)) is True
    assert _ask_recursive(Q.diagonal(V1.T*(X + Z)*V1)) is True
    assert _ask_recursive(Q.diagonal(MatrixSlice(Y, (0, 1), (1, 2)))) is True
    assert _ask_recursive(Q.diagonal(V1.T*(V1 + V2))) is True
    assert _ask_recursive(Q.diagonal(X**3), Q.diagonal(X)) is True
    assert _ask_recursive(Q.diagonal(Identity(3))) is True
    assert _ask_recursive(Q.diagonal(DiagMatrix(V1))) is True
    assert _ask_recursive(Q.diagonal(DiagonalMatrix(X))) is True


def test_non_atoms():
    assert _ask_recursive(Q.real(Trace(X)), Q.positive(Trace(X))) is True

@XFAIL
def test_non_trivial_implies():
    X = MatrixSymbol('X', 3, 3)
    Y = MatrixSymbol('Y', 3, 3)
    assert _ask_recursive(Q.lower_triangular(X+Y), Q.lower_triangular(X) &
               Q.lower_triangular(Y)) is True
    assert _ask_recursive(Q.triangular(X), Q.lower_triangular(X)) is True
    assert _ask_recursive(Q.triangular(X+Y), Q.lower_triangular(X) &
               Q.lower_triangular(Y)) is True

def test_MatrixSlice():
    X = MatrixSymbol('X', 4, 4)
    B = MatrixSlice(X, (1, 3), (1, 3))
    C = MatrixSlice(X, (0, 3), (1, 3))
    assert _ask_recursive(Q.symmetric(B), Q.symmetric(X)) is True
    assert _ask_recursive(Q.invertible(B), Q.invertible(X)) is True
    assert _ask_recursive(Q.diagonal(B), Q.diagonal(X)) is True
    assert _ask_recursive(Q.orthogonal(B), Q.orthogonal(X)) is True
    assert _ask_recursive(Q.upper_triangular(B), Q.upper_triangular(X)) is True

    assert _ask_recursive(Q.symmetric(C), Q.symmetric(X)) is None
    assert _ask_recursive(Q.invertible(C), Q.invertible(X)) is None
    assert _ask_recursive(Q.diagonal(C), Q.diagonal(X)) is None
    assert _ask_recursive(Q.orthogonal(C), Q.orthogonal(X)) is None
    assert _ask_recursive(Q.upper_triangular(C), Q.upper_triangular(X)) is None

def test_det_trace_positive():
    X = MatrixSymbol('X', 4, 4)
    assert _ask_recursive(Q.positive(Trace(X)), Q.positive_definite(X)) is True
    assert _ask_recursive(Q.positive(Determinant(X)), Q.positive_definite(X)) is True

def test_field_assumptions():
    X = MatrixSymbol('X', 4, 4)
    Y = MatrixSymbol('Y', 4, 4)
    assert _ask_recursive(Q.real_elements(X), Q.real_elements(X)) is True
    assert _ask_recursive(Q.integer_elements(X), Q.real_elements(X)) is None
    assert _ask_recursive(Q.complex_elements(X), Q.real_elements(X)) is True
    assert _ask_recursive(Q.complex_elements(X**2), Q.real_elements(X)) is True
    assert _ask_recursive(Q.real_elements(X**2), Q.integer_elements(X)) is True
    assert _ask_recursive(Q.real_elements(X+Y), Q.real_elements(X)) is None
    assert _ask_recursive(Q.real_elements(X+Y), Q.real_elements(X) & Q.real_elements(Y)) is True
    from sympy.matrices.expressions.hadamard import HadamardProduct
    assert _ask_recursive(Q.real_elements(HadamardProduct(X, Y)),
                    Q.real_elements(X) & Q.real_elements(Y)) is True
    assert _ask_recursive(Q.complex_elements(X+Y), Q.real_elements(X) & Q.complex_elements(Y)) is True

    assert _ask_recursive(Q.real_elements(X.T), Q.real_elements(X)) is True
    assert _ask_recursive(Q.real_elements(X.I), Q.real_elements(X) & Q.invertible(X)) is True
    assert _ask_recursive(Q.real_elements(Trace(X)), Q.real_elements(X)) is True
    assert _ask_recursive(Q.integer_elements(Determinant(X)), Q.integer_elements(X)) is True
    assert _ask_recursive(Q.integer_elements(X.I), Q.integer_elements(X)) is None
    alpha = Symbol('alpha')
    assert _ask_recursive(Q.real_elements(alpha*X), Q.real_elements(X) & Q.real(alpha)) is True
    assert _ask_recursive(Q.real_elements(LofLU(X)), Q.real_elements(X)) is True
    e = Symbol('e', integer=True, negative=True)
    assert _ask_recursive(Q.real_elements(X**e), Q.real_elements(X) & Q.invertible(X)) is True
    assert _ask_recursive(Q.real_elements(X**e), Q.real_elements(X)) is None

def test_matrix_element_sets():
    X = MatrixSymbol('X', 4, 4)
    assert _ask_recursive(Q.real(X[1, 2]), Q.real_elements(X)) is True
    assert _ask_recursive(Q.integer(X[1, 2]), Q.integer_elements(X)) is True
    assert _ask_recursive(Q.complex(X[1, 2]), Q.complex_elements(X)) is True
    assert _ask_recursive(Q.integer_elements(Identity(3))) is True
    assert _ask_recursive(Q.integer_elements(ZeroMatrix(3, 3))) is True
    assert _ask_recursive(Q.integer_elements(OneMatrix(3, 3))) is True
    from sympy.matrices.expressions.fourier import DFT
    assert _ask_recursive(Q.complex_elements(DFT(3))) is True


def test_matrix_element_sets_slices_blocks():
    X = MatrixSymbol('X', 4, 4)
    assert _ask_recursive(Q.integer_elements(X[:, 3]), Q.integer_elements(X)) is True
    assert _ask_recursive(Q.integer_elements(BlockMatrix([[X], [X]])),
                        Q.integer_elements(X)) is True

def test_matrix_element_sets_determinant_trace():
    assert _ask_recursive(Q.integer(Determinant(X)), Q.integer_elements(X)) is True
    assert _ask_recursive(Q.integer(Trace(X)), Q.integer_elements(X)) is True
