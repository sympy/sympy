from sympy import (
    Field, Set, AdditionOperator, MultiplicationOperator
)

A = Set('A')
a, b, e1, e2 = [A.element(n) for n in ('a', 'b', 'e1', 'e2')]
add = AdditionOperator(A**2, A, e1)
mul = MultiplicationOperator(A**2, A, e2)
F = Field('F', (A,), (add, mul))

def test_Ring():
    assert F.div(a, b) == F.mul(a, F.pow(b, -1))
    assert F.div(a, e2, evaluate=True) == a
    assert F.div(e2, a, evaluate=True) == F.pow(a, -1)
    assert F.div(a, a, evaluate=True) == e2
