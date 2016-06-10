from sympy.matrices import Matrix
from sympy.matrices.expressions.dotproduct import DotProduct
from sympy.utilities.pytest import raises

A = Matrix(3, 1, [1, 2, 3])
B = Matrix(3, 1, [1, 3, 5])
C = Matrix(4, 1, [1, 2, 4, 5])


def test_docproduct():
    assert DotProduct(A, B).doit() == 22
    raises(TypeError, lambda: DotProduct(A.T, B).doit())
    raises(TypeError, lambda: DotProduct(A, B.T).doit())
    raises(TypeError, lambda: DotProduct(B, C).doit())
    assert DotProduct(A.T, B.T).doit() == 22
