from sympy import Set, Map, BinaryOperator
from sympy.algebras import (
    Magma, Semigroup, Quasigroup, Monoid,
)
from sympy.testing.pytest import raises

A = Set('A')
B = Set('B', (A,))

f1 = Map('f2', domain=A, codomain=A)
f2 = Map('f1', domain=A**2, codomain=A)

class F3_1(BinaryOperator):
    is_left_divisible=True
    domain = A**2
    codomain = A
f3_1 = F3_1()

class F3_2(BinaryOperator):
    is_right_divisible=True
    domain = A**2
    codomain = A
f3_2 = F3_2()

class F4(BinaryOperator):
    is_associative=True
    domain = A**2
    codomain = A
f4 = F4()

def test_Magma():
    # magma consists of one set and one binary operator
    raises(TypeError, lambda: Magma('M', (A,B), (f2,)))
    raises(TypeError, lambda: Magma('M', (A,), (f1,)))
    raises(TypeError, lambda: Magma('M', (A,), (f2,f2)))

def test_Semigroup():
    # Semigroup's operator must be associative
    raises(TypeError, lambda: Semigroup('SG', (A,), (f2,)))

def test_Quasigroup():
    # Quasigroup's operator must have left and right division
    raises(TypeError, lambda: Quasigroup('Q', (A,), (f3_1,)))
    raises(TypeError, lambda: Quasigroup('Q', (A,), (f3_2,)))

def test_Monoid():
    # Monoid's operator must have identity
    raises(TypeError, lambda: Monoid('M', (A,), (f4,)))
