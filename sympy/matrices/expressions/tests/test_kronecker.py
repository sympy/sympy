from __future__ import annotations
from sympy.core.mod import Mod
from sympy.core.numbers import I
from sympy.core.symbol import symbols
from sympy.functions.elementary.integers import floor
from sympy.matrices.dense import (Matrix, eye)
from sympy.matrices.immutable import ImmutableMatrix
from sympy.matrices import MatrixSymbol, Identity
from sympy.matrices.expressions import det, trace

from sympy.matrices.expressions.kronecker import (KroneckerProduct,
                                                  kronecker_product,
                                                  combine_kronecker)


mat1 = Matrix([[1, 2 * I], [1 + I, 3]])
mat2 = Matrix([[2 * I, 3], [4 * I, 2]])

i, j, k, n, m, o, p, x = symbols('i,j,k,n,m,o,p,x')
Z = MatrixSymbol('Z', n, n)
W = MatrixSymbol('W', m, m)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', n, m)
C = MatrixSymbol('C', m, k)


def test_KroneckerProduct():
    assert isinstance(KroneckerProduct(A, B), KroneckerProduct)
    assert KroneckerProduct(A, B).subs(A, C) == KroneckerProduct(C, B)
    assert KroneckerProduct(A, C).shape == (n*m, m*k)
    assert (KroneckerProduct(A, C) + KroneckerProduct(-A, C)).is_ZeroMatrix
    assert (KroneckerProduct(W, Z) * KroneckerProduct(W.I, Z.I)).is_Identity


def test_KroneckerProduct_identity():
    assert KroneckerProduct(Identity(m), Identity(n)) == Identity(m*n)
    assert KroneckerProduct(eye(2), eye(3)) == eye(6)


def test_KroneckerProduct_explicit():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    kp = KroneckerProduct(X, Y)
    assert kp.shape == (4, 4)
    assert kp.as_explicit() == Matrix(
        [
            [X[0, 0]*Y[0, 0], X[0, 0]*Y[0, 1], X[0, 1]*Y[0, 0], X[0, 1]*Y[0, 1]],
            [X[0, 0]*Y[1, 0], X[0, 0]*Y[1, 1], X[0, 1]*Y[1, 0], X[0, 1]*Y[1, 1]],
            [X[1, 0]*Y[0, 0], X[1, 0]*Y[0, 1], X[1, 1]*Y[0, 0], X[1, 1]*Y[0, 1]],
            [X[1, 0]*Y[1, 0], X[1, 0]*Y[1, 1], X[1, 1]*Y[1, 0], X[1, 1]*Y[1, 1]]
        ]
    )


def test_tensor_product_adjoint():
    assert KroneckerProduct(I*A, B).adjoint() == \
        -I*KroneckerProduct(A.adjoint(), B.adjoint())
    assert KroneckerProduct(mat1, mat2).adjoint() == \
        kronecker_product(mat1.adjoint(), mat2.adjoint())


def test_tensor_product_conjugate():
    assert KroneckerProduct(I*A, B).conjugate() == \
        -I*KroneckerProduct(A.conjugate(), B.conjugate())
    assert KroneckerProduct(mat1, mat2).conjugate() == \
        kronecker_product(mat1.conjugate(), mat2.conjugate())


def test_tensor_product_transpose():
    assert KroneckerProduct(I*A, B).transpose() == \
        I*KroneckerProduct(A.transpose(), B.transpose())
    assert KroneckerProduct(mat1, mat2).transpose() == \
        kronecker_product(mat1.transpose(), mat2.transpose())


def test_KroneckerProduct_is_associative():
    assert kronecker_product(A, kronecker_product(
        B, C)) == kronecker_product(kronecker_product(A, B), C)
    assert kronecker_product(A, kronecker_product(
        B, C)) == KroneckerProduct(A, B, C)


def test_KroneckerProduct_is_bilinear():
    assert kronecker_product(x*A, B) == x*kronecker_product(A, B)
    assert kronecker_product(A, x*B) == x*kronecker_product(A, B)


def test_KroneckerProduct_determinant():
    kp = kronecker_product(W, Z)
    assert det(kp) == det(W)**n * det(Z)**m


def test_KroneckerProduct_trace():
    from sympy.matrices.expressions.trace import Trace
    from sympy.matrices.exceptions import NonSquareMatrixError
    from sympy.testing.pytest import raises

    # Square factors evaluate to product of individual traces
    kp = kronecker_product(W, Z)
    assert trace(kp) == trace(W)*trace(Z)

    # Rectangular factors producing square KroneckerProduct
    C_rect = MatrixSymbol('C_rect', 2, 3)
    D_rect = MatrixSymbol('D_rect', 3, 2)
    kp_sq = KroneckerProduct(C_rect, D_rect)
    assert kp_sq.shape == (6, 6)
    assert kp_sq.is_square is True
    assert isinstance(trace(kp_sq), Trace)
    assert Trace(kp_sq).doit() == Trace(kp_sq)

    # Numerical equivalence with explicit matrix trace
    A_num = Matrix([[1, 2, 3], [4, 5, 6]])
    B_num = Matrix([[7, 8], [9, 10], [11, 12]])
    explicit_val = trace(kronecker_product(A_num, B_num))
    assert Trace(kp_sq).as_explicit().subs({C_rect: A_num, D_rect: B_num}).doit() == explicit_val

    # Reversed rectangular factors
    kp_rev = KroneckerProduct(D_rect, C_rect)
    assert Trace(kp_rev).as_explicit().subs({C_rect: A_num, D_rect: B_num}).doit() == trace(kronecker_product(B_num, A_num))

    # Multi-factor KroneckerProduct with rectangular arguments
    Y_sq = MatrixSymbol('Y_sq', 2, 2)
    kp_multi = KroneckerProduct(Y_sq, C_rect, D_rect)
    assert kp_multi.shape == (12, 12)
    assert kp_multi.is_square is True
    assert isinstance(trace(kp_multi), Trace)

    # Multi-factor all square
    kp_all_sq = KroneckerProduct(Y_sq, W, Z)
    assert trace(kp_all_sq) == trace(Y_sq) * trace(W) * trace(Z)

    # Known rectangular KroneckerProduct raises NonSquareMatrixError
    kp_non_sq = KroneckerProduct(C_rect, Y_sq)
    assert kp_non_sq.shape == (4, 6)
    assert kp_non_sq.is_square is False
    raises(NonSquareMatrixError, lambda: trace(kp_non_sq))
    raises(NonSquareMatrixError, lambda: Trace(kp_non_sq))

    # Indeterminate symbolic dimensions leave Trace unevaluated
    kp_indet = KroneckerProduct(C_rect, W)
    assert kp_indet.is_square is None
    assert isinstance(trace(kp_indet), Trace)


def test_KroneckerProduct_isnt_commutative():
    assert KroneckerProduct(A, B) != KroneckerProduct(B, A)
    assert KroneckerProduct(A, B).is_commutative is False


def test_KroneckerProduct_extracts_commutative_part():
    assert kronecker_product(x * A, 2 * B) == x * \
        2 * KroneckerProduct(A, B)


def test_KroneckerProduct_inverse():
    kp = kronecker_product(W, Z)
    assert kp.inverse() == kronecker_product(W.inverse(), Z.inverse())


def test_KroneckerProduct_combine_add():
    kp1 = kronecker_product(A, B)
    kp2 = kronecker_product(C, W)
    assert combine_kronecker(kp1*kp2) == kronecker_product(A*C, B*W)


def test_KroneckerProduct_combine_mul():
    X = MatrixSymbol('X', m, n)
    Y = MatrixSymbol('Y', m, n)
    # A sum of Kronecker products is a single Kronecker product only when
    # the factors coincide in all slots except at most one; merging every
    # slot at once used to produce the wrong result (A + B) x (X + Y):
    kp1 = kronecker_product(A, X)
    kp2 = kronecker_product(B, Y)
    assert combine_kronecker(kp1 + kp2) == kp1 + kp2

    # One differing slot merges by linearity:
    assert combine_kronecker(kronecker_product(A, X) + kronecker_product(B, X)) == \
        kronecker_product(A + B, X)
    assert combine_kronecker(kronecker_product(A, X) + kronecker_product(A, Y)) == \
        kronecker_product(A, X + Y)

    # Scalar multiples of Kronecker products are merged as well:
    assert combine_kronecker(2*kronecker_product(A, X) + 3*kronecker_product(B, X)) == \
        kronecker_product(2*A + 3*B, X)
    assert combine_kronecker(2*kronecker_product(A, X) + 3*kronecker_product(A, X)) == \
        5*kronecker_product(A, X)

    # More than two addends:
    assert combine_kronecker(kronecker_product(A, X) + kronecker_product(B, X)
                             + kronecker_product(A, Y)) in [
        kronecker_product(A + B, X) + kronecker_product(A, Y),
        kronecker_product(A, X + Y) + kronecker_product(B, X)]


def test_KroneckerProduct_combine_add_numeric():
    # Numeric check that combine_kronecker preserves the value of sums:
    m1 = Matrix([[1, 2], [3, 4]])
    m2 = Matrix([[0, 1], [1, 0]])
    m3 = Matrix([[2, 0], [0, 2]])
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    W2 = MatrixSymbol('W2', 2, 2)
    expr = KroneckerProduct(X, Y) + KroneckerProduct(W2, Y)
    combined = combine_kronecker(expr)
    reps = {X: m1, Y: m2, W2: m3}
    assert expr.subs(reps).doit().as_explicit() == \
        combined.subs(reps).doit().as_explicit()


def test_KroneckerProduct_combine_pow():
    X = MatrixSymbol('X', n, n)
    Y = MatrixSymbol('Y', n, n)
    assert combine_kronecker(KroneckerProduct(
        X, Y)**x) == KroneckerProduct(X**x, Y**x)
    assert combine_kronecker(x * KroneckerProduct(X, Y)
                             ** 2) == x * KroneckerProduct(X**2, Y**2)
    assert combine_kronecker(
        x * (KroneckerProduct(X, Y)**2) * KroneckerProduct(A, B)) == x * KroneckerProduct(X**2 * A, Y**2 * B)
    # cannot simplify because of non-square arguments to kronecker product:
    assert combine_kronecker(KroneckerProduct(A, B.T) ** m) == KroneckerProduct(A, B.T) ** m


def test_KroneckerProduct_expand():
    X = MatrixSymbol('X', n, n)
    Y = MatrixSymbol('Y', n, n)

    assert KroneckerProduct(X + Y, Y + Z).expand(kroneckerproduct=True) == \
        KroneckerProduct(X, Y) + KroneckerProduct(X, Z) + \
        KroneckerProduct(Y, Y) + KroneckerProduct(Y, Z)

def test_KroneckerProduct_entry():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', o, p)

    assert KroneckerProduct(A, B)._entry(i, j) == A[Mod(floor(i/o), n), Mod(floor(j/p), m)]*B[Mod(i, o), Mod(j, p)]


def test_KroneckerProduct_singletons():
    A = MatrixSymbol('A', 3, 3)
    A1 = MatrixSymbol('A1', 1, 1)
    B = MatrixSymbol('B', 3, 3)
    B1 = MatrixSymbol('B1', 1, 1)

    assert KroneckerProduct(A, A1, B, B1).doit() == KroneckerProduct(A, B, A1*B1)
    assert KroneckerProduct(A1, A, B1).doit() == KroneckerProduct(A, A1*B1)
    assert KroneckerProduct(A1, B1).doit() == A1*B1
    assert KroneckerProduct(A1, A).doit() == KroneckerProduct(A1, A)
    assert KroneckerProduct(Identity(1), A, A1).doit() == KroneckerProduct(A, A1)
    assert KroneckerProduct(2*A1, A, 3*B1).doit() == 6*KroneckerProduct(A, A1*B1)

    expr = KroneckerProduct(A, A1, B, B1)
    assert expr.doit().as_explicit() == expr.as_explicit()

    p, q = symbols('p, q', commutative=False)
    P = ImmutableMatrix([[p]])
    Q = ImmutableMatrix([[q]])
    assert KroneckerProduct(A, P, B, Q).doit() == KroneckerProduct(A, P, B, Q)
    assert KroneckerProduct(A, P, B, A1).doit() == KroneckerProduct(A, P, B, A1)
    assert KroneckerProduct(P, A1, A, B1).doit() == KroneckerProduct(P, A, A1*B1)


def test_KroneckerProduct_zeromatrix():
    from sympy.matrices.expressions.special import ZeroMatrix
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, k)
    Z_sym = ZeroMatrix(n, m)
    Z_sq = ZeroMatrix(m, m)

    # Unevaluated KroneckerProduct preserves arguments but knows is_zero_matrix
    kp1 = KroneckerProduct(A, Z_sq)
    assert kp1.is_zero_matrix is True
    assert kp1.is_ZeroMatrix is False
    # Generic symbolic product: zero status is unknown (None), not False
    assert KroneckerProduct(A, B).is_zero_matrix is None

    # doit() canonicalization evaluates to ZeroMatrix
    assert kp1.doit() == ZeroMatrix(n*m, m**2)
    assert KroneckerProduct(Z_sym, B).doit() == ZeroMatrix(n*m, m*k)
    assert KroneckerProduct(Z_sym, Z_sq).doit() == ZeroMatrix(n*m, m**2)

    # kronecker_product helper function evaluates to ZeroMatrix
    assert kronecker_product(A, Z_sq) == ZeroMatrix(n*m, m**2)

    # Explicit concrete matrices
    M = Matrix([[1, 2], [3, 4]])
    Z_exp = Matrix.zeros(2, 2)
    assert kronecker_product(M, Z_exp) == Matrix.zeros(4, 4)
    assert KroneckerProduct(M, ZeroMatrix(2, 2)).doit() == ZeroMatrix(4, 4)


def test_KroneckerProduct_power():
    from sympy.matrices.expressions.special import Identity, ZeroMatrix
    from sympy.matrices.expressions.matpow import MatPow

    p = symbols('p', integer=True)
    x = symbols('x')

    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', m, m)
    C = MatrixSymbol('C', k, k)
    K = KroneckerProduct(A, B)

    # Square factor power evaluations
    assert K**0 == Identity(n*m)
    assert K**1 == K
    assert K**2 == KroneckerProduct(A**2, B**2)
    assert K**5 == KroneckerProduct(A**5, B**5)
    assert K**(-1) == KroneckerProduct(A**(-1), B**(-1))
    assert K**(-2) == KroneckerProduct(A**(-2), B**(-2))
    assert K**p == KroneckerProduct(A**p, B**p)

    # Multiple factors
    assert KroneckerProduct(A, B, C)**3 == KroneckerProduct(A**3, B**3, C**3)

    # Non-integer powers remain un-evaluated MatPow
    K_pow_x = K**x
    assert isinstance(K_pow_x, MatPow) and K_pow_x.base == K

    # Square KroneckerProduct with rectangular constituent factors does not distribute
    R1 = MatrixSymbol('R1', 2, 4)
    R2 = MatrixSymbol('R2', 2, 1)
    K_rect = KroneckerProduct(R1, R2)
    assert K_rect.is_square is True
    K_rect_sq = K_rect**2
    assert isinstance(K_rect_sq, MatPow) and K_rect_sq.base == K_rect

    # Explicit concrete matrices
    M1 = Matrix([[1, 2], [3, 4]])
    M2 = Matrix([[2, 0], [0, 1]])
    K_conc = KroneckerProduct(M1, M2)
    assert (K_conc**2).as_explicit() == (K_conc.as_explicit())**2
    assert (K_conc**(-2)).as_explicit() == (K_conc.as_explicit())**(-2)

    # Identity and ZeroMatrix combinations
    assert KroneckerProduct(Identity(n), Identity(m))**2 == Identity(n*m)
    assert KroneckerProduct(ZeroMatrix(n, n), A)**2 == ZeroMatrix(n**2, n**2)


def test_kronecker_combine_add_with_nonkrons():
    # Issue #30900 (item E2-8): combine_kronecker failed with TypeError
    # when adding KroneckerProduct and non-Kronecker terms
    F = MatrixSymbol('F', n*m, n*m)
    kp1 = KroneckerProduct(A, B)
    kp2 = KroneckerProduct(B, B)
    assert combine_kronecker(kp1 + kp2 + F) == KroneckerProduct(A + B, B) + F
    assert combine_kronecker(kp1 + F) == kp1 + F


def test_kronecker_extract_commutative_matmul():
    # Issue #30900 (item E2-6): KroneckerProduct(2*A*B, C).doit() failed
    # because extract_commutative built scalar Mul on non-commutative matrix factors
    D = MatrixSymbol('D', k, n)
    kp1 = KroneckerProduct(2*A*C, B)
    assert kp1.doit() == 2 * KroneckerProduct(A*C, B)

    kp2 = KroneckerProduct(2*A*C, 3*C*D)
    assert kp2.doit() == 6 * KroneckerProduct(A*C, C*D)

    kp3 = KroneckerProduct(2*A, 3*B)
    assert kp3.doit() == 6 * KroneckerProduct(A, B)


def test_kronecker_combine_mul_partial_matching():
    # Chained Kronecker products where initial pairs don't match shape,
    # but subsequent adjacent pairs do
    X1 = MatrixSymbol('X1', 2, 2)
    X2 = MatrixSymbol('X2', 2, 2)
    X3 = MatrixSymbol('X3', 2, 2)
    A1 = MatrixSymbol('A1', 2, 3)
    A2 = MatrixSymbol('A2', 4, 5)
    B1 = MatrixSymbol('B1', 3, 2)
    B2 = MatrixSymbol('B2', 5, 4)

    kp_X = KroneckerProduct(X1, X2, X3)
    kp_A = KroneckerProduct(A1, A2)
    kp_B = KroneckerProduct(B1, B2)

    res = combine_kronecker(kp_X * kp_A * kp_B)
    assert res == kp_X * KroneckerProduct(A1*B1, A2*B2)
