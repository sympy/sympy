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
    associative = True
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
    assert M.sub(r, s) == M.add(r, add.inverse_element(s))
    assert M.sub(x, y) == M.add(x, op.inverse_element(y))
    assert M.sub(x, x, evaluate=True) == e
    assert M.add(M.mul(r, x), M.mul(s, x), evaluate=True) == M.mul(M.add(r, s), x)
    assert M.mul(r, s, y) == M.mul(M.mul(r, s), y)
    assert M.mul(r, r, y, evaluate=True) == M.mul(M.pow(r,2), y)
