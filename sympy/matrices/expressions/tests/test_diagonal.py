from __future__ import annotations
from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices.expressions.diagonal import DiagonalMatrix, DiagonalOf, DiagMatrix, diagonalize_vector
from sympy.assumptions.ask import (Q, ask)
from sympy.core.symbol import Symbol
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.matrices.dense import Matrix
from sympy.core.singleton import S
from sympy.matrices.expressions.matmul import MatMul
from sympy.matrices.expressions.special import Identity, ZeroMatrix
from sympy.matrices.expressions.trace import Trace, trace
from sympy.matrices.expressions.determinant import Determinant
from sympy.matrices.expressions.inverse import Inverse
from sympy.matrices.expressions.hadamard import HadamardPower
from sympy.matrices.expressions.matpow import MatPow
from sympy.matrices.expressions.conjugate import Conjugate
from sympy.matrices.expressions.adjoint import Adjoint
from sympy.matrices.expressions.transpose import Transpose
from sympy.matrices.exceptions import NonSquareMatrixError
from sympy.testing.pytest import raises


n = Symbol('n')
m = Symbol('m')


def test_DiagonalMatrix():
    x = MatrixSymbol('x', n, m)
    D = DiagonalMatrix(x)
    assert D.diagonal_length is None
    assert D.shape == (n, m)

    x = MatrixSymbol('x', n, n)
    D = DiagonalMatrix(x)
    assert D.diagonal_length == n
    assert D.shape == (n, n)
    assert D[1, 2] == 0
    assert D[1, 1] == x[1, 1]
    i = Symbol('i')
    j = Symbol('j')
    x = MatrixSymbol('x', 3, 3)
    ij = DiagonalMatrix(x)[i, j]
    assert ij != 0
    assert ij.subs({i:0, j:0}) == x[0, 0]
    assert ij.subs({i:0, j:1}) == 0
    assert ij.subs({i:1, j:1}) == x[1, 1]
    assert ask(Q.diagonal(D))  # affirm that D is diagonal

    x = MatrixSymbol('x', n, 3)
    D = DiagonalMatrix(x)
    assert D.diagonal_length == 3
    assert D.shape == (n, 3)
    assert D[2, m] == KroneckerDelta(2, m)*x[2, m]
    assert D[3, m] == 0
    raises(IndexError, lambda: D[m, 3])

    x = MatrixSymbol('x', 3, n)
    D = DiagonalMatrix(x)
    assert D.diagonal_length == 3
    assert D.shape == (3, n)
    assert D[m, 2] == KroneckerDelta(m, 2)*x[m, 2]
    assert D[m, 3] == 0
    raises(IndexError, lambda: D[3, m])

    x = MatrixSymbol('x', n, m)
    D = DiagonalMatrix(x)
    assert D.diagonal_length is None
    assert D.shape == (n, m)
    assert D[m, 4] != 0

    x = MatrixSymbol('x', 3, 4)
    assert [DiagonalMatrix(x)[i] for i in range(12)] == [
        x[0, 0], 0, 0, 0, 0, x[1, 1], 0, 0, 0, 0, x[2, 2], 0]

    # shape is retained, issue 12427
    assert (
        DiagonalMatrix(MatrixSymbol('x', 3, 4))*
        DiagonalMatrix(MatrixSymbol('x', 4, 2))).shape == (3, 2)


def test_DiagonalOf():
    x = MatrixSymbol('x', n, n)
    d = DiagonalOf(x)
    assert d.shape == (n, 1)
    assert d.diagonal_length == n
    assert d[2, 0] == d[2] == x[2, 2]

    x = MatrixSymbol('x', n, m)
    d = DiagonalOf(x)
    assert d.shape == (None, 1)
    assert d.diagonal_length is None
    assert d[2, 0] == d[2] == x[2, 2]

    d = DiagonalOf(MatrixSymbol('x', 4, 3))
    assert d.shape == (3, 1)
    d = DiagonalOf(MatrixSymbol('x', n, 3))
    assert d.shape == (3, 1)
    d = DiagonalOf(MatrixSymbol('x', 3, n))
    assert d.shape == (3, 1)
    x = MatrixSymbol('x', n, m)
    assert [DiagonalOf(x)[i] for i in range(4)] ==[
        x[0, 0], x[1, 1], x[2, 2], x[3, 3]]


def test_DiagMatrix():
    x = MatrixSymbol('x', n, 1)
    d = DiagMatrix(x)
    assert d.shape == (n, n)
    assert d[0, 1] == 0
    assert d[0, 0] == x[0, 0]

    a = MatrixSymbol('a', 1, 1)
    d = diagonalize_vector(a)
    assert isinstance(d, MatrixSymbol)
    assert a == d
    assert diagonalize_vector(Identity(3)) == Identity(3)
    assert DiagMatrix(Identity(3)).doit() == Identity(3)
    assert isinstance(DiagMatrix(Identity(3)), DiagMatrix)

    # A diagonal matrix is equal to its transpose:
    assert DiagMatrix(x).T == DiagMatrix(x)
    assert diagonalize_vector(x.T) == DiagMatrix(x)

    dx = DiagMatrix(x)
    assert dx[0, 0] == x[0, 0]
    assert dx[1, 1] == x[1, 0]
    assert dx[0, 1] == 0
    assert dx[0, m] == x[0, 0]*KroneckerDelta(0, m)

    z = MatrixSymbol('z', 1, n)
    dz = DiagMatrix(z)
    assert dz[0, 0] == z[0, 0]
    assert dz[1, 1] == z[0, 1]
    assert dz[0, 1] == 0
    assert dz[0, m] == z[0, m]*KroneckerDelta(0, m)

    v = MatrixSymbol('v', 3, 1)
    dv = DiagMatrix(v)
    assert dv.as_explicit() == Matrix([
        [v[0, 0], 0, 0],
        [0, v[1, 0], 0],
        [0, 0, v[2, 0]],
    ])

    v = MatrixSymbol('v', 1, 3)
    dv = DiagMatrix(v)
    assert dv.as_explicit() == Matrix([
        [v[0, 0], 0, 0],
        [0, v[0, 1], 0],
        [0, 0, v[0, 2]],
    ])

    dv = DiagMatrix(3*v)
    assert dv.args == (3*v,)
    assert dv.doit() == 3*DiagMatrix(v)
    assert isinstance(dv.doit(), MatMul)

    a = MatrixSymbol("a", 3, 1).as_explicit()
    expr = DiagMatrix(a)
    result = Matrix([
        [a[0, 0], 0, 0],
        [0, a[1, 0], 0],
        [0, 0, a[2, 0]],
    ])
    assert expr.doit() == result
    expr = DiagMatrix(a.T)
    assert expr.doit() == result


def test_DiagMatrix_zeromatrix():
    # Issue #30900 (item E3-7): DiagMatrix(ZeroMatrix(3, 1)).doit() mutated
    # the square matrix shape (3, 3) into a non-square shape (3, 1).
    z_col = ZeroMatrix(3, 1)
    d_col = DiagMatrix(z_col)
    assert d_col.shape == (3, 3)
    assert d_col.doit() == ZeroMatrix(3, 3)
    assert d_col.doit().shape == (3, 3)
    assert diagonalize_vector(z_col) == ZeroMatrix(3, 3)

    z_row = ZeroMatrix(1, 3)
    d_row = DiagMatrix(z_row)
    assert d_row.shape == (3, 3)
    assert d_row.doit() == ZeroMatrix(3, 3)
    assert d_row.doit().shape == (3, 3)
    assert diagonalize_vector(z_row) == ZeroMatrix(3, 3)

    # Transposed zero vector
    d_trans = DiagMatrix(z_col.T)
    assert d_trans.shape == (3, 3)
    assert d_trans.doit() == ZeroMatrix(3, 3)

    # Symbolic dimensions
    k = Symbol('k', integer=True)
    z_sym = ZeroMatrix(k, 1)
    d_sym = DiagMatrix(z_sym)
    assert d_sym.shape == (k, k)
    assert d_sym.doit() == ZeroMatrix(k, k)
    assert d_sym.doit().shape == (k, k)
    assert diagonalize_vector(z_sym) == ZeroMatrix(k, k)

    z_sym_row = ZeroMatrix(1, k)
    d_sym_row = DiagMatrix(z_sym_row)
    assert d_sym_row.shape == (k, k)
    assert d_sym_row.doit() == ZeroMatrix(k, k)
    assert d_sym_row.doit().shape == (k, k)


def test_DiagMatrix_calculus():
    v = MatrixSymbol('v', 3, 1)
    d = DiagMatrix(v)

    # Trace
    assert Trace(d).doit() == v[0, 0] + v[1, 0] + v[2, 0]
    v_row = MatrixSymbol('v_row', 1, 3)
    assert Trace(DiagMatrix(v_row)).doit() == v_row[0, 0] + v_row[0, 1] + v_row[0, 2]
    assert Trace(DiagMatrix(ZeroMatrix(0, 1))).doit() == S.Zero
    # Symbolic dimension leaves Trace unevaluated
    k = Symbol('k')
    v_sym = MatrixSymbol('v_sym', k, 1)
    assert Trace(DiagMatrix(v_sym)).doit() == Trace(DiagMatrix(v_sym))

    # Determinant
    assert Determinant(d).doit() == v[0, 0] * v[1, 0] * v[2, 0]
    assert Determinant(DiagMatrix(v_row)).doit() == v_row[0, 0] * v_row[0, 1] * v_row[0, 2]
    assert Determinant(DiagMatrix(ZeroMatrix(0, 1))).doit() == S.One
    assert Determinant(DiagMatrix(v_sym)).doit() == Determinant(DiagMatrix(v_sym))

    # Power & Inverse
    assert d**0 == Identity(3)
    assert d**1 == d
    assert (d**2).doit() == DiagMatrix(HadamardPower(v, 2))
    assert (d**-1).doit() == DiagMatrix(HadamardPower(v, -1))
    assert (d**-2).doit() == DiagMatrix(HadamardPower(v, -2))
    assert Inverse(d).doit() == DiagMatrix(HadamardPower(v, -1))
    # Power composition
    assert ((d**2)**3).doit() == DiagMatrix(HadamardPower(v, 6))
    # Symbolic power stays unevaluated MatPow
    assert (d**k).doit() == MatPow(d, k)
    # Empty vector doit
    assert DiagMatrix(Matrix(0, 1, [])).doit() == Matrix(0, 0, [])

    # Conjugate & Adjoint
    assert d.conjugate() == DiagMatrix(Conjugate(v))
    assert d.adjoint() == DiagMatrix(Conjugate(v))


def test_DiagonalMatrix_calculus():
    M = MatrixSymbol('M', 3, 3)
    dm = DiagonalMatrix(M)

    # Trace
    assert Trace(dm).doit() == Trace(M)
    assert trace(dm) == Trace(M)
    # Rectangular trace
    M_rect = MatrixSymbol('M_rect', 2, 3)
    dm_rect = DiagonalMatrix(M_rect)
    raises(NonSquareMatrixError, lambda: Trace(dm_rect))

    # Determinant
    assert Determinant(dm).doit() == M[0, 0] * M[1, 1] * M[2, 2]
    M_0 = MatrixSymbol('M_0', 0, 0)
    assert Determinant(DiagonalMatrix(M_0)).doit() == S.One
    # Symbolic dimension
    k = Symbol('k')
    M_sym = MatrixSymbol('M_sym', k, k)
    assert Determinant(DiagonalMatrix(M_sym)).doit() == Determinant(DiagonalMatrix(M_sym))

    # Transpose, Conjugate, Adjoint
    assert dm.T == DiagonalMatrix(Transpose(M))
    assert dm.conjugate() == DiagonalMatrix(Conjugate(M))
    assert dm.adjoint() == DiagonalMatrix(Adjoint(M))

    # Power & Inverse
    assert dm**0 == Identity(3)
    assert dm**1 == dm
    assert dm**2 == DiagMatrix(Matrix([M[0, 0]**2, M[1, 1]**2, M[2, 2]**2]))
    assert dm**-1 == DiagMatrix(Matrix([M[0, 0]**-1, M[1, 1]**-1, M[2, 2]**-1]))
    assert dm**-2 == DiagMatrix(Matrix([M[0, 0]**-2, M[1, 1]**-2, M[2, 2]**-2]))
    assert (dm**2).doit() == Matrix([
        [M[0, 0]**2, 0, 0],
        [0, M[1, 1]**2, 0],
        [0, 0, M[2, 2]**2],
    ])
    assert (dm**-2).doit() == Matrix([
        [M[0, 0]**-2, 0, 0],
        [0, M[1, 1]**-2, 0],
        [0, 0, M[2, 2]**-2],
    ])
    assert Inverse(dm).doit() == DiagMatrix(Matrix([
        [M[0, 0]**-1],
        [M[1, 1]**-1],
        [M[2, 2]**-1],
    ]))
    assert Inverse(dm).doit().doit() == Matrix([
        [M[0, 0]**-1, 0, 0],
        [0, M[1, 1]**-1, 0],
        [0, 0, M[2, 2]**-1],
    ])
    assert (dm**k).doit() == MatPow(dm, k)
    # ZeroMatrix and Identity powers
    assert (DiagonalMatrix(ZeroMatrix(3, 3))**2).doit() == ZeroMatrix(3, 3)
    assert (DiagonalMatrix(Identity(3))**5).doit() == Identity(3)

    # as_explicit on concrete Matrix
    M_dense = Matrix([[1, 2], [3, 4]])
    dm_dense = DiagonalMatrix(M_dense)
    assert dm_dense.as_explicit() == Matrix([[1, 0], [0, 4]])


def test_diagonal_assumptions():
    v = MatrixSymbol('v', 3, 1)
    d = DiagMatrix(v)
    M = MatrixSymbol('M', 3, 3)
    dm = DiagonalMatrix(M)
    M_rect = MatrixSymbol('M_rect', 2, 3)
    dm_rect = DiagonalMatrix(M_rect)

    # Symmetric
    assert ask(Q.symmetric(d)) is True
    assert ask(Q.symmetric(dm)) is True
    assert ask(Q.symmetric(dm_rect)) is False

    # Triangular
    assert ask(Q.lower_triangular(d)) is True
    assert ask(Q.upper_triangular(d)) is True
    assert ask(Q.lower_triangular(dm)) is True
    assert ask(Q.upper_triangular(dm)) is True

    # Element predicates
    assert ask(Q.integer_elements(d), Q.integer_elements(v)) is True
    assert ask(Q.integer_elements(dm), Q.integer_elements(M)) is True

    assert ask(Q.real_elements(d), Q.real_elements(v)) is True
    assert ask(Q.real_elements(dm), Q.real_elements(M)) is True

    assert ask(Q.complex_elements(d), Q.complex_elements(v)) is True
    assert ask(Q.complex_elements(dm), Q.complex_elements(M)) is True
