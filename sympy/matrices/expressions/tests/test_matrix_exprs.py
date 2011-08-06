from sympy.utilities.pytest import raises
from sympy import S, symbols, Symbol
from sympy.matrices import (eye, MatrixSymbol, Transpose, Inverse, ShapeError,
        MatMul, Identity, BlockMatrix, BlockDiagMatrix, block_collapse, Matrix,
        ZeroMatrix, MatAdd)

def test_transpose():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)

    assert Transpose(A).shape == (m,n)
    assert Transpose(A*B).shape == (l,n)
    assert Transpose(Transpose(A)) == A

    assert Transpose(eye(3)) == eye(3)

def test_inverse():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', n, n)
    D = MatrixSymbol('D', n, n)

    raises(ShapeError, "Inverse(A)")
    assert Inverse(Inverse(C)) == C

    assert Inverse(C)*C == Identity(C.n)

    assert Inverse(eye(3)) == eye(3)

    assert Inverse(S(3)) == S(1)/3

    assert Inverse(Identity(n)) == Identity(n)

    assert Inverse(C*D) == D.I*C.I

def test_shape():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    assert A.shape == (n, m)
    assert (A*B).shape == (n, l)
    raises(ShapeError, 'B*A')

def test_matexpr():
    n, m, l = symbols('n m l', integer=True)
    x = Symbol('x')
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)

    assert (x*A).shape == A.shape
    assert (x*A).__class__ == MatMul
    assert 2*A - A - A == ZeroMatrix(*A.shape)

def test_subs():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', m, l)

    assert A.subs(n,m).shape == (m,m)

    assert (A*B).subs(B,C) == A*C

    assert (A*B).subs(l,n).is_square


def test_BlockMatrix():
    n,m,l,k = symbols('n m l k', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, k)
    C = MatrixSymbol('C', l, m)
    D = MatrixSymbol('D', l, k)
    X = BlockMatrix(Matrix([[A,B],[C,D]]))


    assert X.shape == (l+n, k+m)
    assert Transpose(X) == BlockMatrix(Matrix([[A.T, C.T], [B.T, D.T]]))
    assert Transpose(X).shape == X.shape[::-1]
    assert X.blockshape == (2,2)

    E = MatrixSymbol('E', m, 1)
    F = MatrixSymbol('F', k, 1)

    Y = BlockMatrix(Matrix([[E], [F]]))

    assert (X*Y).shape == (l+n, 1)
    assert block_collapse(X*Y)[0,0] == A*E + B*F
    assert block_collapse(X*Y)[1,0] == C*E + D*F
    assert Transpose(block_collapse(Transpose(X*Y))) == block_collapse(X*Y)

def test_squareBlockMatrix():
    n,m,l,k = symbols('n m l k', integer=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, m)
    C = MatrixSymbol('C', m, n)
    D = MatrixSymbol('D', m, m)
    X = BlockMatrix([[A,B],[C,D]])
    Y = BlockMatrix([[A]])


    assert X.is_square

    assert block_collapse(X+Identity(m+n)) == BlockMatrix(
        [[A+Identity(n), B], [C, D+Identity(m)]])

    assert (X + MatrixSymbol('Q', n+m, n+m)).is_Add
    assert (X * MatrixSymbol('Q', n+m, n+m)).is_Mul

    assert Y.I[0,0] == A.I
    assert Inverse(X, expand=True) == BlockMatrix([
        [(-B*D.I*C + A).I, -A.I*B*(D+-C*A.I*B).I],
        [-(D-C*A.I*B).I*C*A.I, (D-C*A.I*B).I]])

    assert Inverse(X, expand=False).__class__ == Inverse


def test_BlockDiagMatrix():
    n,m,l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', m, m)
    C = MatrixSymbol('C', l, l)

    X = BlockDiagMatrix(A,B,C)

    assert X[1,1] == B
    assert X.shape == (n+m+l, n+m+l)
    assert all(X[i,j].is_ZeroMatrix if i!=j else X[i,j] in [A,B,C]
            for i in range(3) for j in range(3))

    assert block_collapse(X.I * X).is_Identity

    assert block_collapse(X*X) == BlockDiagMatrix(A**2, B**2, C**2)

    assert block_collapse(X+X) == BlockDiagMatrix(2*A, 2*B, 2*C)

def test_ZeroMatrix():
    n,m = symbols('n m', integer=True)
    A = MatrixSymbol('A', n, m)
    Z = ZeroMatrix(n, m)

    assert A+Z == A
    assert A*Z.T == ZeroMatrix(n,n)
    assert Z*A.T == ZeroMatrix(n,n)
    assert A-A == ZeroMatrix(*A.shape)

    assert Transpose(Z) == ZeroMatrix(m, n)

def test_Identity():
    n,m = symbols('n m', integer=True)
    A = MatrixSymbol('A', n, m)
    In = Identity(n)
    Im = Identity(m)

    assert A*Im == A
    assert In*A == A

    assert Transpose(In) == In
    assert Inverse(In) == In


def test_MatAdd():
    n, m = symbols('n m', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, m)


    assert (A+B).shape == A.shape
    assert MatAdd(A, -A, 2*B).is_Mul

    raises(ShapeError, "A + B.T")
    raises(ValueError, "A+1")

def test_MatMul():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)

    assert (2*A*B).shape == (n,l)

    raises(ShapeError, "B*A")
    assert (2*A).shape == A.shape

    assert MatMul(A, ZeroMatrix(m,m), B) == ZeroMatrix(n,l)

def test_MatPow():
    n = Symbol('n', integer=True)
    A = MatrixSymbol('A', n, n)

    assert Inverse(A).is_Pow
    assert (A*A).is_Pow
    assert (A*A).exp == 2
    assert (A*A).base == A
    assert (A**n).exp == n

    assert A**0 == Identity(n)
    assert A**1 == A
    raises(ShapeError, "MatrixSymbol('B', 3,2)**2")

def test_linear_factors():
    from sympy.matrices import MatrixSymbol, linear_factors
    n, m, l = symbols('n m l')
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', n, l)
    assert linear_factors(2*A*B + C, B, C) == { C: Identity(n), B: 2*A}
