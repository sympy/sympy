from sympy.matrices.expressions import MatrixSymbol, HadamardProduct
from sympy.matrices import ShapeError
from sympy import symbols
from sympy.utilities.pytest import raises


def test_HadamardProduct():
    n, m = symbols('n,m')
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, m)

    assert HadamardProduct(A, B, A).shape == A.shape

    raises(ShapeError, lambda: HadamardProduct(A, B.T))
    raises(TypeError,  lambda: A + 1)
    raises(TypeError,  lambda: 5 + A)
    raises(TypeError,  lambda: 5 - A)

    assert HadamardProduct(A, 2*B, -A)[1, 1] == -2 * A[1, 1]**2 * B[1, 1]
