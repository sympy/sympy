from sympy import Set, BinaryOperator
from sympy.algebras import Monoid
from sympy.testing.pytest import raises

A = Set('A')

class F(BinaryOperator):
    is_associative=True
    domain = A**2
    codomain = A
f = F()

def test_Monoid():
    # Monoid's operator must have identity
    raises(TypeError, lambda: Monoid('M', (A,), (f,)))
