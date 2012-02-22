from __future__ import division
from sympy.linalg.dokmatrix import DOKMatrix

def test_getitem():
    A = DOKMatrix(((1,2,0),(4,0,6),(7,8,0)))
    assert A[0, 0] == 1
    assert A[0, 2] == 0

def test_setitem():
    A = DOKMatrix(((1,2,3),(0,0,0),(3,4,0)))
    A[0, 0] = 10
    A[1, 1] = 15
    assert A == DOKMatrix(((10,2,3),(0,15,0),(3,4,0)))

def test_copyinmatrix():
    A = DOKMatrix(((1,0,3),(0,5,6),(7,0,9)))
    B = DOKMatrix(((10,11),(12,13)))
    keys = (slice(None, 2),slice(None, 2))
    A.copyin_matrix(keys, B)
    assert A == DOKMatrix(((10,11,3),(12,13,6),(7,0,9)))

def test_transpose():
    A = DOKMatrix(((1,0,3),(0,5,6),(7,0,9)))
    assert A.T == DOKMatrix(((1,0,7),(0,5,0),(3,6,9)))

def test_add():
    A = DOKMatrix(((1,0,3),(0,5,6),(7,0,9)))
    B = DOKMatrix(((1,2,3),(0,0,0),(3,4,0)))
    assert A + B == DOKMatrix(((2,2,6),(0,5,6),(10,4,9)))

def test_mul():
    A = DOKMatrix(((1,2,3),(4,0,1),(9,8,0)))
    B = DOKMatrix(((0,0,1),(2,1,0),(0,2,7)))
    C = DOKMatrix(((4,8,22),(0,2,11),(16,8,9)))
    assert A * B == C

def test_submatrix():
    A = DOKMatrix(((1,2,3),(4,0,1),(9,8,0)))
    B = DOKMatrix(((0,0,1),(2,1,0),(0,2,7)))
    C = DOKMatrix(((4,8,27),(0,2,11),(16,8,9)))
    keys = (slice(None, 2), slice(None, 2))
    A.submatrix(keys) == DOKMatrix(((1,2),(4,0)))
    B.submatrix(keys) == DOKMatrix(((0,0),(2,1)))
    C.submatrix(keys) == DOKMatrix(((4,8),(0,2)))

def test_solve():
    A = DOKMatrix(((1,2,3),(4,0,1),(9,8,0)))
    rhs = DOKMatrix([[3,4,5]]).T
    for method in ["LDL", "CH"]:
        assert A * A.solve(rhs, method=method) == rhs, method
