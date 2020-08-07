from sympy import (
Set, AdditionOperator, BinaryOperator, Ring,
VectorAdditionOperator, AbelianGroup, ScalarMultiplicationOperator,
Module
)

A = Set('A')
r, s, e1, e2 = [A.element(n) for n in ('r', 's', 'e1', 'e2')]
add = AdditionOperator(A**2, A, e1)
class MonoidOp(BinaryOperator):
    name = '*'
    domain = A*A
    codomain = A
    is_associative = True
    identity = e2
mul = MonoidOp()
R = Ring('R', (A,), (add, mul))

X = Set('X')
x, y, e = [X.element(i) for i in 'xye']
op = VectorAdditionOperator(X**2, X, e)
G = AbelianGroup('G', (X,), (op,))

smul = ScalarMultiplicationOperator(R*G, G)

M = Module('M', (R, G), (smul,))

def test_Module():
    assert M.add(r, r, evaluate=True) == M.mul(M.add(e2, e2), r)
    assert M.add(x, x, evaluate=True) == M.mul(M.add(e2, e2), x)
