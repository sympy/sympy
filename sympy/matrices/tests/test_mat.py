from sympy.matrices.mat import doktocsr
from sympy import SparseMatrix
a = SparseMatrix(3, 4, {(0, 0): 1, (0, 1): 2, (1, 1): 3, (1, 2): 9, (2, 1): 1, (2, 2): 4})
b = SparseMatrix(4, 6, {(0, 0): 10, (0, 1): 20, (1, 1): 30, (1, 3): 40, (2, 2): 50, (2, 3): 60, (2, 4): 70, (3, 5): 80})

def  test_doktocsr():
    assert doktocsr(a) == [[1, 2, 3, 9, 1, 4], [0, 1, 1, 2, 1, 2], [0, 2, 4, 6]]
    assert doktocsr(b) == [[10, 20, 30, 40, 50, 60, 70, 80], [0, 1, 1, 3, 2, 3, 4, 5], [0, 2, 4, 7, 8]]
