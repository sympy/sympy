from sympy import S
from sympy.sets.undefinedset import UndefinedSet

A = UndefinedSet('A')
B = UndefinedSet('B', (A,))
C = UndefinedSet('C', (B,))

def test_UndefinedSet():
    assert A == Set('A')
    assert B == Set('B', (A,))

    assert B.is_subset(A)
    assert A.is_superset(B)
    assert C.is_subset(A)
    assert A.is_superset(C)

    assert A.intersection(B) == B
    assert B.intersection(A) == B
    assert A.union(B) == A
    assert B.union(A) == A
