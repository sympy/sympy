from sympy import Set, Map
from sympy.algebras import Magma
from sympy.testing.pytest import raises

A = Set('A')
B = Set('B', (A,))

f1 = Map('f1', domain=A, codomain=A)
f2 = Map('f2', domain=A**2, codomain=A)

def test_Magma():
    # magma consists of one set and one binary operator
    raises(TypeError, lambda: Magma('M', (A,B), (f2,)))
    raises(TypeError, lambda: Magma('M', (A,), (f1,)))
    raises(TypeError, lambda: Magma('M', (A,), (f2,f2)))
