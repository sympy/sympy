from sympy import Set, BinaryOperator
from sympy.algebras import Loop
from sympy.testing.pytest import raises

A = Set('A')

class F(BinaryOperator):
    is_left_divisible = is_right_divisible = True
    domain = A**2
    codomain = A
f = F()

def test_Loop():
    # Loop's operator must have identity
    raises(TypeError, lambda: Loop('L', (A,), (f,)))
