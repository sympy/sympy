from sympy import Set, Map
from sympy.algebras import AlgebraicStructure

A = Set('A')
B = Set('B', (A,))
f_A = Map('f', domain=A, codomain=A)
f_B = Map('f', domain=B, codomain=B)

S_A = AlgebraicStructure('S_A', (A,), (f_A,))
S_B = AlgebraicStructure('S_B', (B,), (f_B,))

a, b = A.element('a'), B.element('b')

def test_AlgebraicStructure():

    # check substructure
    assert S_A.is_substructure(S_A)
    assert S_A.is_superstructure(S_A)
    assert S_B.is_substructure(S_A)
    assert S_A.is_superstructure(S_B)

    # check element
    assert a in S_A
    assert a not in S_B
    assert b in S_A
    assert b in S_B

    # check codomain
    assert f_A(a) in S_A
    assert f_B(b) in S_B
    assert f_A(a) not in S_B
    assert f_B(b) in S_A
