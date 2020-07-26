from sympy.core.expr import unchanged
from sympy.core.symbol import symbols
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.matrices.expressions.matadd import MatAdd
from sympy.matrices.expressions.special import (
    ZeroMatrix, GenericZeroMatrix, Identity, GenericIdentity)
from sympy.matrices.expressions.matmul import MatMul
from sympy.testing.pytest import raises


def test_zero_matrix_creation():
    assert unchanged(ZeroMatrix, 2, 2)
    assert unchanged(ZeroMatrix, 0, 0)
    raises(ValueError, lambda: ZeroMatrix(-1, 2))
    raises(ValueError, lambda: ZeroMatrix(2.0, 2))
    raises(ValueError, lambda: ZeroMatrix(2j, 2))
    raises(ValueError, lambda: ZeroMatrix(2, -1))
    raises(ValueError, lambda: ZeroMatrix(2, 2.0))
    raises(ValueError, lambda: ZeroMatrix(2, 2j))

    n = symbols('n')
    assert unchanged(ZeroMatrix, n, n)
    n = symbols('n', integer=False)
    raises(ValueError, lambda: ZeroMatrix(n, n))
    n = symbols('n', negative=True)
    raises(ValueError, lambda: ZeroMatrix(n, n))


def test_generic_zero_matrix():
    z = GenericZeroMatrix()
    n = symbols('n', integer=True)
    A = MatrixSymbol("A", n, n)

    assert z == z
    assert z != A
    assert A != z

    assert z.is_ZeroMatrix

    raises(TypeError, lambda: z.shape)
    raises(TypeError, lambda: z.rows)
    raises(TypeError, lambda: z.cols)

    assert MatAdd() == z
    assert MatAdd(z, A) == MatAdd(A)
    # Make sure it is hashable
    hash(z)


def test_identity_matrix_creation():
    assert Identity(2)
    assert Identity(0)
    raises(ValueError, lambda: Identity(-1))
    raises(ValueError, lambda: Identity(2.0))
    raises(ValueError, lambda: Identity(2j))

    n = symbols('n')
    assert Identity(n)
    n = symbols('n', integer=False)
    raises(ValueError, lambda: Identity(n))
    n = symbols('n', negative=True)
    raises(ValueError, lambda: Identity(n))


def test_generic_identity():
    I = GenericIdentity()
    n = symbols('n', integer=True)
    A = MatrixSymbol("A", n, n)

    assert I == I
    assert I != A
    assert A != I

    assert I.is_Identity
    assert I**-1 == I

    raises(TypeError, lambda: I.shape)
    raises(TypeError, lambda: I.rows)
    raises(TypeError, lambda: I.cols)

    assert MatMul() == I
    assert MatMul(I, A) == MatMul(A)
    # Make sure it is hashable
    hash(I)
