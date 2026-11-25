from __future__ import annotations
from itertools import product

from sympy.core.relational import (Equality, Unequality)
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.integrals.integrals import integrate
from sympy.matrices.dense import (Matrix, eye, zeros)
from sympy.matrices.immutable import ImmutableMatrix
from sympy.matrices import SparseMatrix
from sympy.matrices.immutable import \
    ImmutableDenseMatrix, ImmutableSparseMatrix
from sympy.abc import x, y
from sympy.testing.pytest import raises
from sympy.core.symbol import symbols
from sympy.core.numbers import Rational
from sympy.polys.domains import QQ

IM = ImmutableDenseMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
ISM = ImmutableSparseMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
ieye = ImmutableDenseMatrix(eye(3))


def test_creation():
    assert IM.shape == ISM.shape == (3, 3)
    assert IM[1, 2] == ISM[1, 2] == 6
    assert IM[2, 2] == ISM[2, 2] == 9


def test_immutability():
    with raises(TypeError):
        IM[2, 2] = 5
    with raises(TypeError):
        ISM[2, 2] = 5


def test_slicing():
    assert IM[1, :] == ImmutableDenseMatrix([[4, 5, 6]])
    assert IM[:2, :2] == ImmutableDenseMatrix([[1, 2], [4, 5]])
    assert ISM[1, :] == ImmutableSparseMatrix([[4, 5, 6]])
    assert ISM[:2, :2] == ImmutableSparseMatrix([[1, 2], [4, 5]])


def test_subs():
    A = ImmutableMatrix([[1, 2], [3, 4]])
    B = ImmutableMatrix([[1, 2], [x, 4]])
    C = ImmutableMatrix([[-x, x*y], [-(x + y), y**2]])
    assert B.subs(x, 3) == A
    assert (x*B).subs(x, 3) == 3*A
    assert (x*eye(2) + B).subs(x, 3) == 3*eye(2) + A
    assert C.subs([[x, -1], [y, -2]]) == A
    assert C.subs([(x, -1), (y, -2)]) == A
    assert C.subs({x: -1, y: -2}) == A
    assert C.subs({x: y - 1, y: x - 1}, simultaneous=True) == \
        ImmutableMatrix([[1 - y, (x - 1)*(y - 1)], [2 - x - y, (x - 1)**2]])


def test_as_immutable():
    data = [[1, 2], [3, 4]]
    X = Matrix(data)
    assert sympify(X) == X.as_immutable() == ImmutableMatrix(data)

    data = {(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}
    X = SparseMatrix(2, 2, data)
    assert sympify(X) == X.as_immutable() == ImmutableSparseMatrix(2, 2, data)


def test_function_return_types():
    # Lets ensure that decompositions of immutable matrices remain immutable
    # I.e. do MatrixBase methods return the correct class?
    X = ImmutableMatrix([[1, 2], [3, 4]])
    Y = ImmutableMatrix([[1], [0]])
    q, r = X.QRdecomposition()
    assert (type(q), type(r)) == (ImmutableMatrix, ImmutableMatrix)

    assert type(X.LUsolve(Y)) == ImmutableMatrix
    assert type(X.QRsolve(Y)) == ImmutableMatrix

    X = ImmutableMatrix([[5, 2], [2, 7]])
    assert X.T == X
    assert X.is_symmetric
    assert type(X.cholesky()) == ImmutableMatrix
    L, D = X.LDLdecomposition()
    assert (type(L), type(D)) == (ImmutableMatrix, ImmutableMatrix)

    X = ImmutableMatrix([[1, 2], [2, 1]])
    assert X.is_diagonalizable()
    assert X.det() == -3
    assert X.norm(2) == 3

    assert type(X.eigenvects()[0][2][0]) == ImmutableMatrix

    assert type(zeros(3, 3).as_immutable().nullspace()[0]) == ImmutableMatrix

    X = ImmutableMatrix([[1, 0], [2, 1]])
    assert type(X.lower_triangular_solve(Y)) == ImmutableMatrix
    assert type(X.T.upper_triangular_solve(Y)) == ImmutableMatrix

    assert type(X.minor_submatrix(0, 0)) == ImmutableMatrix

# issue 6279
# https://github.com/sympy/sympy/issues/6279
# Test that Immutable _op_ Immutable => Immutable and not MatExpr


def test_immutable_evaluation():
    X = ImmutableMatrix(eye(3))
    A = ImmutableMatrix(3, 3, range(9))
    assert isinstance(X + A, ImmutableMatrix)
    assert isinstance(X * A, ImmutableMatrix)
    assert isinstance(X * 2, ImmutableMatrix)
    assert isinstance(2 * X, ImmutableMatrix)
    assert isinstance(A**2, ImmutableMatrix)


def test_deterimant():
    assert ImmutableMatrix(4, 4, lambda i, j: i + j).det() == 0


def test_Equality():
    assert Equality(IM, IM) is S.true
    assert Unequality(IM, IM) is S.false
    assert Equality(IM, IM.subs(1, 2)) is S.false
    assert Unequality(IM, IM.subs(1, 2)) is S.true
    assert Equality(IM, 2) is S.false
    assert Unequality(IM, 2) is S.true
    M = ImmutableMatrix([x, y])
    assert Equality(M, IM) is S.false
    assert Unequality(M, IM) is S.true
    assert Equality(M, M.subs(x, 2)).subs(x, 2) is S.true
    assert Unequality(M, M.subs(x, 2)).subs(x, 2) is S.false
    assert Equality(M, M.subs(x, 2)).subs(x, 3) is S.false
    assert Unequality(M, M.subs(x, 2)).subs(x, 3) is S.true


def test_integrate():
    intIM = integrate(IM, x)
    assert intIM.shape == IM.shape
    assert all(intIM[i, j] == (1 + j + 3*i)*x for i, j in
                product(range(3), range(3)))



def test_hash():
    A = ImmutableMatrix([[1, 0]])
    B = ImmutableSparseMatrix([[1, 0]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# 1x1 matrices

def test_hash_1x1_numeric():
    A = ImmutableMatrix([[1]])
    B = ImmutableSparseMatrix([[1]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_1x1_symbolic():
    x = symbols('x')

    A = ImmutableMatrix([[x]])
    B = ImmutableSparseMatrix([[x]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_1x1_mixed():
    x = symbols('x')

    A = ImmutableMatrix([[x + 1]])
    B = ImmutableSparseMatrix([[x + 1]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# 1x2 matrices

def test_hash_1x2_numeric():
    A = ImmutableMatrix([[1, 2]])
    B = ImmutableSparseMatrix([[1, 2]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_1x2_symbolic():
    x, y = symbols('x y')

    A = ImmutableMatrix([[x, y]])
    B = ImmutableSparseMatrix([[x, y]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_1x2_mixed():
    x = symbols('x')

    A = ImmutableMatrix([[x, 1]])
    B = ImmutableSparseMatrix([[x, 1]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# 2x1 matrices

def test_hash_2x1_numeric():
    A = ImmutableMatrix([[1], [2]])
    B = ImmutableSparseMatrix([[1], [2]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_2x1_symbolic():
    x, y = symbols('x y')

    A = ImmutableMatrix([[x], [y]])
    B = ImmutableSparseMatrix([[x], [y]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_2x1_mixed():
    x = symbols('x')

    A = ImmutableMatrix([[x], [1]])
    B = ImmutableSparseMatrix([[x], [1]])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# 3x3 matrices

def test_hash_3x3_numeric():
    A = ImmutableMatrix([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ])

    B = ImmutableSparseMatrix([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_3x3_symbolic():
    x, y, z = symbols('x y z')

    A = ImmutableMatrix([
        [x, y, z],
        [y, z, x],
        [z, x, y]
    ])

    B = ImmutableSparseMatrix([
        [x, y, z],
        [y, z, x],
        [z, x, y]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_3x3_mixed():
    x, y, z = symbols('x y z')

    A = ImmutableMatrix([
        [x, 1, 0],
        [2, y, 3],
        [0, 4, z]
    ])

    B = ImmutableSparseMatrix([
        [x, 1, 0],
        [2, y, 3],
        [0, 4, z]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# 3x2 matrices

def test_hash_3x2_numeric():
    A = ImmutableMatrix([
        [1, 2],
        [3, 4],
        [5, 6]
    ])

    B = ImmutableSparseMatrix([
        [1, 2],
        [3, 4],
        [5, 6]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_3x2_symbolic():
    x, y, z = symbols('x y z')

    A = ImmutableMatrix([
        [x, y],
        [y, z],
        [z, x]
    ])

    B = ImmutableSparseMatrix([
        [x, y],
        [y, z],
        [z, x]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_3x2_mixed():
    x, y, z = symbols('x y z')

    A = ImmutableMatrix([
        [x, 1],
        [2, y],
        [z, 0]
    ])

    B = ImmutableSparseMatrix([
        [x, 1],
        [2, y],
        [z, 0]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Identity matrix

def test_hash_identity():
    A = ImmutableMatrix.eye(3)
    B = ImmutableSparseMatrix.eye(3)

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Diagonal matrices

def test_hash_diagonal_symbolic():
    x, y, z = symbols('x y z')

    A = ImmutableMatrix.diag(x, y, z)
    B = ImmutableSparseMatrix.diag(x, y, z)

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_diagonal_mixed():
    x, y = symbols('x y')

    A = ImmutableMatrix.diag(x, 1, y)
    B = ImmutableSparseMatrix.diag(x, 1, y)

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Zero matrix

def test_hash_zero():
    A = ImmutableMatrix.zeros(3, 3)
    B = ImmutableSparseMatrix.zeros(3, 3)

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Empty matrices

def test_hash_empty_0x0():
    A = ImmutableMatrix([])
    B = ImmutableSparseMatrix([])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_empty_0x3():
    A = ImmutableMatrix(0, 3, [])
    B = ImmutableSparseMatrix(0, 3, [])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


def test_hash_empty_3x0():
    A = ImmutableMatrix(3, 0, [])
    B = ImmutableSparseMatrix(3, 0, [])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Negative numbers

def test_hash_negative():
    A = ImmutableMatrix([
        [-1, -2],
        [-3, -4]
    ])

    B = ImmutableSparseMatrix([
        [-1, -2],
        [-3, -4]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Rational numbers

def test_hash_rational():
    A = ImmutableMatrix([
        [Rational(1, 2), 0],
        [0, Rational(3, 4)]
    ])

    B = ImmutableSparseMatrix([
        [Rational(1, 2), 0],
        [0, Rational(3, 4)]
    ])

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Different insertion order

def test_hash_insertion_order():
    A = ImmutableSparseMatrix(2, 2, {
        (0, 0): 1,
        (1, 1): 2
    })

    B = ImmutableSparseMatrix(2, 2, {
        (1, 1): 2,
        (0, 0): 1
    })

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Different domains: Dense vs Sparse

def test_hash_different_domains_dense_sparse():
    A = ImmutableMatrix([[1, 0]])
    B = ImmutableSparseMatrix._fromrep(
        A._rep.convert_to(QQ)
    )

    assert A._rep.domain != B._rep.domain
    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Different domains: Dense vs Dense

def test_hash_different_domains_dense_dense():
    A = ImmutableMatrix([[1, 0]])
    B = ImmutableMatrix._fromrep(
        A._rep.convert_to(QQ)
    )

    assert A._rep.domain != B._rep.domain
    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1


# Different domains: Sparse vs Sparse

def test_hash_different_domains_sparse_sparse():
    A = ImmutableSparseMatrix([[1, 0]])
    B = ImmutableSparseMatrix._fromrep(
        A._rep.convert_to(QQ)
    )

    assert A._rep.domain != B._rep.domain
    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1

def test_hash_different_domains_rational():
    from sympy.polys.domains import EXRAW

    A = ImmutableMatrix([[Rational(1, 2), 0]])
    B = ImmutableSparseMatrix._fromrep(
        A._rep.convert_to(EXRAW)
    )

    assert A == B
    assert hash(A) == hash(B)
    assert len({A, B}) == 1
