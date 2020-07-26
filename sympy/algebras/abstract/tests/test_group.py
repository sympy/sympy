from sympy import Set, BinaryOperator
from sympy.algebras import Group, AbelianGroup
from sympy.testing.pytest import raises

A = Set('A')
e = A.element('e')

class F1(BinaryOperator):
    is_left_divisible = is_right_divisible = True
    identity = e
    domain = A**2
    codomain = A
f1 = F1()

class F2(BinaryOperator):
    is_associative = True
    identity = e
    domain = A**2
    codomain = A
f2 = F2()

class G(BinaryOperator):
    is_associative = True
    is_left_divisible = is_right_divisible = True
    identity = e
    domain = A**2
    codomain = A
g = G()

def test_Group():
    # Group's operator must have left and right division
    raises(TypeError, lambda: Group('G', (A,), (f1,)))
    # Group's operator must be associative
    raises(TypeError, lambda: Group('G', (A,), (f2,)))

def test_AbelianGroup():
    # Abelian group's operator must be commutative
    raises(TypeError, lambda: AbelianGroup('G', (A,), (g,)))
