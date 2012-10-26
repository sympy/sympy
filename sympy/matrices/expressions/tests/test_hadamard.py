from sympy.matrices.expressions import MatrixSymbol, HadamardProduct
from sympy.matrices import ShapeError
from sympy import symbols
from sympy.utilities.pytest import raises


def test_HadamardProduct():
    n, m, k = symbols('n,m,k')
    Z = MatrixSymbol('Z', n, n)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, m)
    C = MatrixSymbol('C', m, k)

    assert HadamardProduct(A, B, A).shape == A.shape

    raises(ShapeError, lambda: HadamardProduct(A, B.T))
    raises(TypeError,  lambda: A + 1)
    raises(TypeError,  lambda: 5 + A)
    raises(TypeError,  lambda: 5 - A)

    assert HadamardProduct(A, 2*B, -A)[1, 1] == -2 * A[1, 1]**2 * B[1, 1]

    mix = HadamardProduct(Z*A, B)*C
    assert mix.shape == (n, k)
