from sympy.matrices.expressions import MatrixSymbol, Adjoint, Trace
from sympy.matrices import eye, Matrix
from sympy import symbols, S, conjugate

n,m,l,k,p = symbols('n m l k p', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)

def test_adjoint():
    Sq = MatrixSymbol('Sq', n, n)

    assert Adjoint(A).shape == (m, n)
    assert Adjoint(A*B).shape == (l, n)
    assert Adjoint(Adjoint(A)) == A

    assert Adjoint(eye(3)) == eye(3)

    assert Adjoint(S(5)) == S(5)

    assert Adjoint(Matrix([[1, 2], [3, 4]])) == Matrix([[1, 3], [2, 4]])

    assert Adjoint(Trace(Sq)) == conjugate(Trace(Sq))
    assert Trace(Adjoint(Sq)) == conjugate(Trace(Sq))

    assert Adjoint(Sq)[0, 1] == conjugate(Sq[1, 0])
