from sympy.matrices.expressions import MatrixSymbol, Adjoint, Trace
from sympy.matrices import eye, Matrix
from sympy import symbols, S, conjugate, adjoint

n,m,l,k,p = symbols('n m l k p', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)

def test_adjoint():
    Sq = MatrixSymbol('Sq', n, n)

    assert Adjoint(A).shape == (m, n)
    assert Adjoint(A*B).shape == (l, n)
    assert adjoint(Adjoint(A)) == A
    assert isinstance(Adjoint(Adjoint(A)), Adjoint)

    assert Adjoint(eye(3)).doit() == eye(3)

    assert Adjoint(S(5)).doit() == S(5)

    assert Adjoint(Matrix([[1, 2], [3, 4]])).doit() == Matrix([[1, 3], [2, 4]])

    assert adjoint(Trace(Sq)) == conjugate(Trace(Sq))
    assert Trace(adjoint(Sq)) == conjugate(Trace(Sq))

    assert Adjoint(Sq)[0, 1] == conjugate(Sq[1, 0])

    assert Adjoint(A*B).doit() == Adjoint(B) * Adjoint(A)
