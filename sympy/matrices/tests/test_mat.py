from sympy.matrices.mat import doktocsr
from sympy import SparseMatrix
a = SparseMatrix(3, 4, {(0, 0): 1, (0, 1): 2, (1, 1): 3, (1, 2): 9, (2, 1): 1, (2, 2): 4})
b = SparseMatrix(4, 6, {(0, 0): 10, (0, 1): 20, (1, 1): 30, (1, 3): 40, (2, 2): 50, (2, 3): 60, (2, 4): 70, (3, 5): 80})
c = SparseMatrix(4, 4, {(1, 1): 12, (1, 3): 2, (2, 0): 15, (2, 2): 12, (3, 3): 4})
d = SparseMatrix(10, 10, {(1, 1): 12, (3, 5): 7, (7, 8): 12})
e = SparseMatrix(3, 3, {(1, 0): 1, (1, 2): 2, (2, 0): 3})

def test_doktocsr():
    assert doktocsr(a) == [[1, 2, 3, 9, 1, 4], [0, 1, 1, 2, 1, 2], [0, 2, 4, 6]]
    assert doktocsr(b) == [[10, 20, 30, 40, 50, 60, 70, 80], [0, 1, 1, 3, 2, 3, 4, 5], [0, 2, 4, 7, 8]]
    assert doktocsr(c) == [[12, 2, 15, 12, 4], [1, 3, 0, 2, 3], [0, 0, 2, 4, 5]]
    assert doktocsr(d) == [[12, 7, 12], [1, 5, 8], [0, 0, 1, 1, 2, 2, 2, 2, 3]]
    assert doktocsr(e) == [[1, 2, 3], [0, 2, 0], [0, 0, 2, 3]]
