from sympy import (
    Ring, Set, BinaryOperator, AdditionOperator
)

A = Set('A')
a, b, e1, e2 = [A.element(n) for n in ('a', 'b', 'e1', 'e2')]
add = AdditionOperator(A**2, A, e1)
class MonoidOp(BinaryOperator):
    domain = A*A
    codomain = A
    is_associative = True
    identity = e2
mul = MonoidOp()
R = Ring('R', (A,), (add, mul))

def test_Ring():
    assert R.add(a, a, evaluate=True) == R.mul(R.add(e2, e2), a)
    assert R.add(a, e1, evaluate=True) == a

    assert R.sub(a, a, evaluate=True) == e1

    assert R.mul(b, e2, evaluate=True) == b
    assert R.pow(b, 1, evaluate=True) == b
    assert R.pow(b, 0, evaluate=True) == e2
