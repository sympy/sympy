from sympy.liealgebras.cartan_matrix import CartanMatrix
from sympy.matrices import Matrix

def test_CartanMatrix_TypeA():
    """Rest of cartan testing should be test_typeX.py"""
    c = CartanMatrix("A3")
    m = Matrix(3, 3, [2, -1, 0, -1, 2, -1, 0, -1, 2])
    assert c == m
