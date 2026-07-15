"""Tests for the .T (transpose) property on all algebraic tensor types."""

from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic import (
    AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor
)


# ---------------------------------------------------------------------------
# AlgebraicZeroTensor transpose tests
# ---------------------------------------------------------------------------

def test_zero_tensor_T_single_factor():
    Z = AlgebraicZeroTensor((3, 4))
    ZT = Z.T
    assert isinstance(ZT, AlgebraicZeroTensor)
    assert ZT.shape == ((4, 3),)


def test_zero_tensor_T_multi_factor():
    Z = AlgebraicZeroTensor(((1, 2), (3, 4)))
    ZT = Z.T
    assert ZT.shape == ((2, 1), (4, 3))


def test_zero_tensor_T_commutativity_pattern():
    Z = AlgebraicZeroTensor(((1, 2), (3, 4)))
    ZT = Z.T
    assert ZT.commutativity_pattern == (1, 1)


def test_zero_tensor_T_double_transpose():
    Z = AlgebraicZeroTensor(((2, 3), (5, 7)))
    assert Z.T.T.shape == Z.shape


# ---------------------------------------------------------------------------
# AlgebraicPureTensor transpose tests
# ---------------------------------------------------------------------------

def test_pure_tensor_T_basic():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    T = AlgebraicPureTensor(A, B)
    TT = T.T
    assert isinstance(TT, AlgebraicPureTensor)
    assert TT.shape == ((4, 3), (5, 4))
    assert TT.factors[0] == A.T
    assert TT.factors[1] == B.T


def test_pure_tensor_T_coeff_preserved():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    T = AlgebraicPureTensor(2, A, B)
    TT = T.T
    assert TT.coeff == 2
    assert TT.shape == ((4, 3), (5, 4))


def test_pure_tensor_T_symbolic_coeff():
    from sympy.abc import x
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    T = AlgebraicPureTensor(x, A, B)
    TT = T.T
    assert TT.coeff == x
    assert TT.shape == ((4, 3), (5, 4))


def test_pure_tensor_T_single_factor():
    A = MatrixSymbol("A", 3, 4)
    # Single factor with coeff=1 unwraps to bare A
    T = AlgebraicPureTensor(A)
    assert T is A
    # Transpose of the unwrapped factor
    T2 = AlgebraicPureTensor(2, A)
    TT = T2.T
    assert TT.coeff == 2
    assert TT.shape == ((4, 3),)


def test_pure_tensor_T_double_transpose():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    T = AlgebraicPureTensor(A, B)
    TT = T.T
    TTT = TT.T
    assert TTT.shape == T.shape


# ---------------------------------------------------------------------------
# AlgebraicTensor transpose tests
# ---------------------------------------------------------------------------

def test_tensor_T_basic():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B),
                        AlgebraicPureTensor(C, D))
    ST = S.T
    assert ST.shape == ((4, 3), (5, 4))


def test_tensor_T_with_zero_anchor():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), Z)
    ST = S.T
    assert ST.shape == ((4, 3), (5, 4))


def test_tensor_T_double_transpose():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B),
                        AlgebraicPureTensor(C, D))
    ST = S.T
    STT = ST.T
    assert STT.shape == S.shape


def test_tensor_T_coefficients_preserved():
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    C = MatrixSymbol("C", 3, 4)
    D = MatrixSymbol("D", 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(2, A, B),
                        AlgebraicPureTensor(3, C, D))
    ST = S.T
    assert ST.shape == ((4, 3), (5, 4))
    for term in ST.terms:
        assert term.coeff in (2, 3)
