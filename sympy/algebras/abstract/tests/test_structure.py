from sympy import Set, Map, S
from sympy.algebras import AlgebraicStructure
from sympy.testing.pytest import raises

A = Set('A')
B = Set('B', (A,))
a, b = A.element('a'), B.element('b')
f_A = Map('f', domain=A**3, codomain=A)
f_B = Map('f', domain=B**3, codomain=B)

S_A = AlgebraicStructure('S_A', (A,), (f_A,))
S_B = AlgebraicStructure('S_B', (B,), (f_B,))

def test_AlgebraicStructure():

    # set must be provided
    raises(TypeError, lambda: AlgebraicStructure('X', (), ()))

    # return set when no operator is given
    assert AlgebraicStructure('X', (S.Reals,), ()) == S.Reals
    assert AlgebraicStructure('X', (S.Integers, S.Reals), ()) == S.Reals

    # check closure
    # domain of operator must be superset of domain of structure
    raises(TypeError, lambda: AlgebraicStructure('X', (A,), (Map('f', domain=B, codomain=A),)))
    AlgebraicStructure('X', (B,), (Map('f', domain=A, codomain=B),)) # not raise error
    # codomain of operator must be subset of domain of structure
    raises(TypeError, lambda: AlgebraicStructure('X', (B,), (Map('f', domain=B, codomain=A),)))
    AlgebraicStructure('X', (A,), (Map('f', domain=A, codomain=B),)) # not raise error

    # element of structure is element of set
    assert S_A.element('a') == A.element('a')

    # element of set is in structure
    assert a in S_A
    assert a not in S_B
    assert b in S_A
    assert b in S_B
    assert f_A(a,a,a) in S_A
    assert f_B(b,b,b) in S_B
    assert f_A(a,a,a) not in S_B
    assert f_B(b,b,b) in S_A

    # check substructure
    assert S_A.is_substructure(S_A)
    assert S_A.is_superstructure(S_A)
    assert S_B.is_substructure(S_A)
    assert S_A.is_superstructure(S_B)

    # check substructure relation with subclasses
    class NewStructure(AlgebraicStructure):
        pass
    ns = NewStructure('NS', (B,), (f_B,))
    assert S_A.is_superstructure(ns)
    assert S_B.is_superstructure(ns)
    assert ns.is_substructure(S_A)
    assert ns.is_substructure(S_B)
