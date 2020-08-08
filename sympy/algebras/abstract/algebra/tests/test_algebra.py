from sympy import (
Set, AdditionOperator, MultiplicationOperator, Field,
VectorAdditionOperator, AbelianGroup, ScalarMultiplicationOperator,
VectorSpace, VectorMultiplicationOperator, Algebra
)

A = Set('A')
a, e1, e2 = [A.element(n) for n in ('a', 'e1', 'e2')]
add = AdditionOperator(A**2, A, e1)
mul = MultiplicationOperator(A**2, A, e2)
F = Field('F', (A,), (add, mul))

X = Set('X')
x, e = [X.element(i) for i in 'xe']
op = VectorAdditionOperator(X**2, X, e)
G = AbelianGroup('G', (X,), (op,))

smul = ScalarMultiplicationOperator(F*G, G)
V = VectorSpace('V', (F, G), (smul,))

vmul = VectorMultiplicationOperator(V*V, G)
A = Algebra('A', (V,), (vmul,))

def test_Algebra():
    assert A.mul(a, x, x, evaluate=True) == smul(a, vmul(x, x))
