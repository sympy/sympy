from sympy import S, Set, Symbol
from sympy.sets.undefinedset import UndefinedSet, SetElement

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

    D = UndefinedSet('D', (S.Reals,))
    assert D.is_subset(S.Reals)
    assert D.is_subset(S.Complexes)
    assert S.Reals.is_superset(D)
    assert S.Complexes.is_superset(D)
    assert D.intersection(S.Reals) == D
    assert S.Reals.intersection(D) == D
    assert D.intersection(S.Complexes) == D
    assert S.Complexes.intersection(D) == D

def test_SetElement():
    c = SetElement('c', C)
    assert c == C.element('c')
    assert (c in C) and (c in B) and (c in A)

    D = Set('D')
    E = Set('E', (A.intersect(D),))
    e = E.element('e')
    assert (e in D) and (e in A)

    assert SetElement('x', S.Complexes) == Symbol('x', complex=True)
    assert SetElement('x', S.Reals) == Symbol('x', real=True)
    assert SetElement('x', S.Rationals) == Symbol('x', rational=True)
    assert SetElement('x', S.Integers) == Symbol('x', integer=True)
    assert SetElement('x', S.Naturals0) == Symbol('x', integer=True, nonnegative=True)
    assert SetElement('x', S.Naturals) == Symbol('x', integer=True, positive=True)
