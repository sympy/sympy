from sympy import Set, Map
from sympy.algebras import (
    Magma, Semigroup
)
from sympy.testing.pytest import raises

A = Set('A')
B = Set('B', (A,))

f = Map('f', domain=A**2, codomain=A)
g = Map('g', domain=A, codomain=A)

def test_Magma():

    # magma consists of one set and one binary operator
    raises(TypeError, lambda: Magma('M', (A,B), (f,)))
    raises(TypeError, lambda: Magma('M', (A,), (g,)))
    raises(TypeError, lambda: Magma('M', (A,), (f,f)))

def test_Semigroup():
    # Semigroup's operator must be associative
    raises(TypeError, lambda: Semigroup('SG', (A,), (f,)))
