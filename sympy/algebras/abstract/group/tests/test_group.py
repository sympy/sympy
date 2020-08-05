from sympy import Set, BinaryOperator
from sympy.algebras import Group, AbelianGroup
from sympy.testing.pytest import raises

A = Set('A')
e = A.element('e')
a, b, e = [A.element(n) for n in 'abe']

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

class G(F1, F2):
    pass
g = G()

class H(G):
    is_commutative = True
h = H()

def test_Group():
    # Group's operator must have left and right division
    raises(TypeError, lambda: Group('G', (A,), (f1,)))
    # Group's operator must be associative
    raises(TypeError, lambda: Group('G', (A,), (f2,)))

    G = Group('G', (A,), (g,))
    G_op = G.operator
    assert G_op(a, G_op(e, b), evaluate=True) == G_op(a, b)
    assert G_op(a, G.inverse(a), evaluate=True) == e

def test_AbelianGroup():
    # Abelian group's operator must be commutative
    raises(TypeError, lambda: AbelianGroup('G', (A,), (g,)))

    G = AbelianGroup('G', (A,), (h,))
    G_op = G.operator
    assert G_op(b, G_op(a, b), evaluate=True) == G_op(a, b, b, evaluate=True)
