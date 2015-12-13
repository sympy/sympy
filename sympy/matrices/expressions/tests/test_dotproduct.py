from sympy.matrices import Matrix
from sympy.matrices.expressions.dotproduct import DotProduct

A = Matrix(3, 1, [1, 2, 3])
B = Matrix(3, 1, [1, 3, 5])

def test_docproduct():
    assert DotProduct(A, B).doit() == 22
