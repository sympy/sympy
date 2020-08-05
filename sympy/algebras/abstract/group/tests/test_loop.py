from sympy import Set, BinaryOperator
from sympy.algebras import Loop
from sympy.testing.pytest import raises

A = Set('A')
a, b, e = [A.element(n) for n in 'abe']

class F(BinaryOperator):
    is_left_divisible = is_right_divisible = True
    domain = A**2
    codomain = A
f = F()

class G(F):
    identity = e
g = G()

def test_Loop():
    # Loop's operator must have identity
    raises(TypeError, lambda: Loop('L', (A,), (f,)))

    L = Loop('L', (A,), (g,))
    L_op = L.operator
    L_ld, L_rd = L.left_division, L.right_division
    assert L_op(a, e, evaluate=True) == a
    assert L_op(e, b, evaluate=True) == b
