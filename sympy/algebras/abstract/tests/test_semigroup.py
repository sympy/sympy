from sympy import Set, Map
from sympy.algebras import Semigroup
from sympy.testing.pytest import raises

A = Set('A')
f = Map('f', domain=A**2, codomain=A)

def test_Semigroup():
    # Semigroup's operator must be associative
    raises(TypeError, lambda: Semigroup('SG', (A,), (f,)))
