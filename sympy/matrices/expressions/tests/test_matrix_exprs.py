from sympy.core import Lambda, S, symbols, Tuple
from sympy.functions import transpose, sin, cos, sqrt
from sympy.simplify import simplify
from sympy.matrices import (
    BlockDiagMatrix, BlockMatrix, block_collapse, eye, FunctionMatrix,
    Identity, ImmutableMatrix, Inverse, MatAdd, MatMul, MatPow, Matrix,
    MatrixExpr, MatrixSymbol, ShapeError, Trace, Transpose, ZeroMatrix,
)
from sympy.utilities.pytest import raises

n, m, l, k, p = symbols('n m l k p', integer=True)
x = symbols('x')
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)
E = MatrixSymbol('E', m, n)


def test_shape():
    assert A.shape == (n, m)
    assert (A*B).shape == (n, l)
    raises(ShapeError, lambda: B*A)


def test_matexpr():
    assert (x*A).shape == A.shape
    assert (x*A).__class__ == MatMul
    assert 2*A - A - A == ZeroMatrix(*A.shape)
    assert (A*B).shape == (n, l)


def test_subs():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', m, l)

    assert A.subs(n, m).shape == (m, m)

    assert (A*B).subs(B, C) == A*C

    assert (A*B).subs(l, n).is_square


def test_BlockMatrix():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, k)
    C = MatrixSymbol('C', l, m)
    D = MatrixSymbol('D', l, k)
    M = MatrixSymbol('M', m + k, p)
    N = MatrixSymbol('N', l + n, k + m)
    X = BlockMatrix(Matrix([[A, B], [C, D]]))

    # block_collapse does nothing on normal inputs
    E = MatrixSymbol('E', n, m)
    assert block_collapse(A + 2*E) == A + 2*E
    F = MatrixSymbol('F', m, m)
    assert block_collapse(E.T*A*F) == E.T*A*F

    assert X.shape == (l + n, k + m)
    assert X.blockshape == (2, 2)
    assert transpose(X) == BlockMatrix(Matrix([[A.T, C.T], [B.T, D.T]]))
    assert transpose(X).shape == X.shape[::-1]

    # Test that BlockMatrices and MatrixSymbols can still mix
    assert (X*M).is_MatMul
    assert X._blockmul(M).is_MatMul
    assert (X*M).shape == (n + l, p)
    assert (X + N).is_MatAdd
    assert X._blockadd(N).is_MatAdd
    assert (X + N).shape == X.shape

    E = MatrixSymbol('E', m, 1)
    F = MatrixSymbol('F', k, 1)

    Y = BlockMatrix(Matrix([[E], [F]]))

    assert (X*Y).shape == (l + n, 1)
    assert block_collapse(X*Y).blocks[0, 0] == A*E + B*F
    assert block_collapse(X*Y).blocks[1, 0] == C*E + D*F

    # block_collapse passes down into container objects, transposes, and inverse
    assert block_collapse(transpose(X*Y)) == transpose(block_collapse(X*Y))
    assert block_collapse(Tuple(X*Y, 2*X)) == (
        block_collapse(X*Y), block_collapse(2*X))

    # Make sure that MatrixSymbols will enter 1x1 BlockMatrix if it simplifies
    Ab = BlockMatrix([[A]])
    Z = MatrixSymbol('Z', *A.shape)
    assert block_collapse(Ab + Z) == A + Z


def test_BlockMatrix_Trace():
    A, B, C, D = map(lambda s: MatrixSymbol(s, 3, 3), 'ABCD')
    X = BlockMatrix([[A, B], [C, D]])
    assert Trace(X) == Trace(A) + Trace(D)


def test_squareBlockMatrix():
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, m)
    C = MatrixSymbol('C', m, n)
    D = MatrixSymbol('D', m, m)
    X = BlockMatrix([[A, B], [C, D]])
    Y = BlockMatrix([[A]])

    assert X.is_square

    assert (block_collapse(X + Identity(m + n)) ==
        BlockMatrix([[A + Identity(n), B], [C, D + Identity(m)]]))
    Q = X + Identity(m + n)

    assert block_collapse(Q.inverse()) == Inverse(block_collapse(Q))

    assert (X + MatrixSymbol('Q', n + m, n + m)).is_MatAdd
    assert (X * MatrixSymbol('Q', n + m, n + m)).is_MatMul

    assert Y.I.blocks[0, 0] == A.I
    assert X.inverse(expand=True) == BlockMatrix([
        [(-B*D.I*C + A).I, -A.I*B*(D + -C*A.I*B).I],
        [-(D - C*A.I*B).I*C*A.I, (D - C*A.I*B).I]])

    assert isinstance(X.inverse(expand=False), Inverse)
    assert isinstance(X.inverse(), Inverse)

    assert not X.is_Identity

    Z = BlockMatrix([[Identity(n), B], [C, D]])
    assert not Z.is_Identity


def test_BlockDiagMatrix():
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', m, m)
    C = MatrixSymbol('C', l, l)
    M = MatrixSymbol('M', n + m + l, n + m + l)

    X = BlockDiagMatrix(A, B, C)
    Y = BlockDiagMatrix(A, 2*B, 3*C)

    assert X.blocks[1, 1] == B
    assert X.shape == (n + m + l, n + m + l)
    assert all(X.blocks[i, j].is_ZeroMatrix if i != j else X.blocks[i, j] in [A, B, C]
            for i in range(3) for j in range(3))

    assert isinstance(block_collapse(X.I * X), Identity)

    assert block_collapse(X*X) == BlockDiagMatrix(A*A, B*B, C*C)
    #XXX: should be == ??
    assert block_collapse(X + X).equals(BlockDiagMatrix(2*A, 2*B, 2*C))
    assert block_collapse(X*Y) == BlockDiagMatrix(A*A, 2*B*B, 3*C*C)
    assert block_collapse(X + Y) == BlockDiagMatrix(2*A, 3*B, 4*C)

    # Ensure that BlockDiagMatrices can still interact with normal MatrixExprs
    assert (X*(2*M)).is_MatMul
    assert (X + (2*M)).is_MatAdd

    assert (X._blockmul(M)).is_MatMul
    assert (X._blockadd(M)).is_MatAdd


def test_ZeroMatrix():
    A = MatrixSymbol('A', n, m)
    Z = ZeroMatrix(n, m)

    assert A + Z == A
    assert A*Z.T == ZeroMatrix(n, n)
    assert Z*A.T == ZeroMatrix(n, n)
    assert A - A == ZeroMatrix(*A.shape)

    assert transpose(Z) == ZeroMatrix(m, n)
    assert Z.conjugate() == Z

    assert ZeroMatrix(n, n)**0 == Identity(n)
    with raises(ShapeError):
        Z**0
    with raises(ShapeError):
        Z**2


def test_Identity():
    A = MatrixSymbol('A', n, m)
    In = Identity(n)
    Im = Identity(m)

    assert A*Im == A
    assert In*A == A

    assert transpose(In) == In
    assert In.inverse() == In
    assert In.conjugate() == In


def test_MatAdd():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, m)

    assert (A + B).shape == A.shape
    assert MatAdd(A, -A, 2*B).is_MatMul

    raises(ShapeError, lambda: A + B.T)
    raises(TypeError, lambda: A + 1)
    raises(TypeError, lambda: 5 + A)
    raises(TypeError, lambda: 5 - A)

    assert MatAdd(A, ZeroMatrix(n, m), -A) == ZeroMatrix(n, m)
    # raises(TypeError, lambda : MatAdd(ZeroMatrix(n,m), S(0)))


def test_multiplication():
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', n, n)

    assert (2*A*B).shape == (n, l)

    assert (A*0*B) == ZeroMatrix(n, l)

    raises(ShapeError, lambda: B*A)
    assert (2*A).shape == A.shape

    assert A * ZeroMatrix(m, m) * B == ZeroMatrix(n, l)

    assert C * Identity(n) * C.I == Identity(n)

    assert B/2 == S.Half*B
    raises(NotImplementedError, lambda: 2/B)

    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, n)
    assert Identity(n) * (A + B) == A + B


def test_MatPow():
    A = MatrixSymbol('A', n, n)

    AA = MatPow(A, 2)
    assert AA.exp == 2
    assert AA.base == A
    assert (A**n).exp == n

    assert A**0 == Identity(n)
    assert A**1 == A
    assert A**2 == AA
    assert A**-1 == Inverse(A)
    assert A**S.Half == sqrt(A)
    raises(ShapeError, lambda: MatrixSymbol('B', 3, 2)**2)


def test_MatrixSymbol():
    n, m, t = symbols('n,m,t')
    X = MatrixSymbol('X', n, m)
    assert X.shape == (n, m)
    raises(TypeError, lambda: MatrixSymbol('X', n, m)(t))  # issue 2756
    assert X.doit() == X


def test_dense_conversion():
    X = MatrixSymbol('X', 2, 2)
    x00, x01, x10, x11 = symbols('X_00 X_01 X_10 X_11')
    assert ImmutableMatrix(X) == ImmutableMatrix([[x00, x01], [x10, x11]])
    assert Matrix(X) == Matrix([[x00, x01], [x10, x11]])


def test_free_symbols():
    assert (C*D).free_symbols == set((C, D))


def test_zero_matmul():
    assert isinstance(S.Zero * MatrixSymbol('X', 2, 2), MatrixExpr)


def test_matadd_simplify():
    A = MatrixSymbol('A', 1, 1)
    assert simplify(MatAdd(A, ImmutableMatrix([[sin(x)**2 + cos(x)**2]]))) == \
        MatAdd(A, ImmutableMatrix([[1]]))


def test_matmul_simplify():
    A = MatrixSymbol('A', 1, 1)
    assert simplify(MatMul(A, ImmutableMatrix([[sin(x)**2 + cos(x)**2]]))) == \
        MatMul(A, ImmutableMatrix([[1]]))
