from sympy import Set, Map
from sympy.algebras import Magma
from sympy.testing.pytest import raises

A = Set('A')
B = Set('B', (A,))

f1 = Map('f1', domain=A, codomain=A)
f2 = Map('f2', domain=A**2, codomain=A)

def test_Magma():
    # magma consists of one set and one binary operator which is closed on the set.
    raises(TypeError, lambda: Magma('M', (A,B), (f2,)))
    raises(TypeError, lambda: Magma('M', (A,), (f1,)))
    raises(TypeError, lambda: Magma('M', (A,), (f2,f2)))

    M = Magma('M', (A,), (f2,))
    M_op = M.operator
    a, b = [A.element(i) for i in 'ab']

    # result of magma operation is closed on the set.
    assert a in M and b in M
    assert M_op(a, b) in M
