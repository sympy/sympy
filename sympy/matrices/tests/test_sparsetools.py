from sympy.matrices.sparsetools import trace, conjugate, transpose, conjugate_transpose
from sympy import SparseMatrix

from sympy import ZZ

def test_trace():
    a = SparseMatrix(7, 8, {(2, 3): ZZ(5), (4, 5):ZZ(12)})
    b = SparseMatrix(10, 10, {(1, 1): ZZ(12), (3, 5): ZZ(7), (7, 8): ZZ(12)})

    assert trace(a, ZZ) == 0
    assert trace(b, ZZ) == 12


def test_conjugate():
    a = SparseMatrix(7, 8, {(2, 3): 5, (4, 5): 12})
    b = SparseMatrix(10, 10, {(1, 1): 12, (3, 5): 7, (7, 8): 12})

    assert conjugate(a, ZZ) == SparseMatrix(7, 8, {(2, 3): 5, (4, 5): 12})
    assert conjugate(b, ZZ) == SparseMatrix(10, 10, {(1, 1): 12, (3, 5): 7, (7, 8): 12})


def test_transpose():
    a = SparseMatrix(7, 8, {(2, 3): ZZ(5), (4, 5):ZZ(12)})
    b = SparseMatrix(10, 10, {(1, 1): ZZ(12), (3, 5): ZZ(7), (7, 8): ZZ(12)})

    assert transpose(a, ZZ) == SparseMatrix(8, 7, {(3, 2): ZZ(5), (5, 4): ZZ(12)})
    assert transpose(b, ZZ) == SparseMatrix(10, 10, {(1, 1): ZZ(12), (5, 3): ZZ(7), (8, 7): ZZ(12)})


def test_conjugate_transpose():
    a = SparseMatrix(7, 8, {(2, 3): 5, (4, 5): 12})
    b = SparseMatrix(10, 10, {(1, 1): 12, (3, 5): 7, (7, 8): 12})

    assert conjugate_transpose(a, ZZ) == SparseMatrix(8, 7, {(3, 2): 5, (5, 4): 12})
    assert conjugate_transpose(b, ZZ) == SparseMatrix(10, 10, {(1, 1): 12, (5, 3): 7, (8, 7): 12})
