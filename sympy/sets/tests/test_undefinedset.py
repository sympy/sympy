from sympy import S
from sympy.sets.undefinedset import UndefinedSet, SetElement

A = UndefinedSet('A')
B = UndefinedSet('B', (A,))
C = UndefinedSet('C', (B,))

def test_UndefinedSet():
    assert B.is_subset(A)
    assert A.is_superset(B)
    assert C.is_subset(A)
    assert A.is_superset(C)
