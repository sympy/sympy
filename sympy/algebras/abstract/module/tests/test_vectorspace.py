from sympy import (
Set, AdditionOperator, MultiplicationOperator, Field,
VectorAdditionOperator, AbelianGroup, ScalarMultiplicationOperator,
VectorSpace
)

A = Set('A')
r, s, e1, e2 = [A.element(n) for n in ('r', 's', 'e1', 'e2')]
add = AdditionOperator(A**2, A, e1)
mul = MultiplicationOperator(A**2, A, e2)
F = Field('F', (A,), (add, mul))

X = Set('X')
x, e = [X.element(i) for i in 'xe']
op = VectorAdditionOperator(X**2, X, e)
G = AbelianGroup('G', (X,), (op,))

smul = ScalarMultiplicationOperator(F*G, G)

V = VectorSpace('V', (F, G), (smul,))

def test_VectorSpace():
    assert V.div(r, s) == V.mul(r, V.pow(s, -1))
    assert V.div(x, s) == V.mul(x, V.pow(s, -1))
    assert V.div(r, e2, evaluate=True) == r
    assert V.div(x, e2, evaluate=True) == x
    assert V.div(e2, r, evaluate=True) == V.pow(r, -1)
