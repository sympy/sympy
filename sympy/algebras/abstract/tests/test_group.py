from sympy import Set, Map
from sympy.algebras import (
    AlgebraicStructure, Magma,
)

A = Set('A')
B = Set('B', (A,))
f_A = Map('f', domain=A*A, codomain=A)
f_B = Map('f', domain=B*B, codomain=B)

M_A = Magma('M_A', (A,), (f_A,))
M_B = Magma('M_B', (B,), (f_B,))

a1, a2 = A.element('a1'), A.element('a2')
b1, b2 = B.element('b1'), B.element('b2')

def test_Magma():

    # check substructure
    assert M_A.is_substructure(M_A)
    assert M_A.is_superstructure(M_A)
    assert M_B.is_substructure(M_A)
    assert M_A.is_superstructure(M_B)

    # check element
    assert a1 in M_A
    assert a1 not in M_B
    assert b1 in M_A
    assert b1 in M_B

    # check codomain
    assert f_A(a1, a2) in M_A
    assert f_B(b1, b2) in M_B
    assert f_A(a1, a2) not in M_B
    assert f_B(b1, b2) in M_A
