from sympy import Set, BinaryOperator
from sympy.algebras import Monoid
from sympy.testing.pytest import raises

A = Set('A')

class F(BinaryOperator):
    is_associative=True
    domain = A**2
    codomain = A
f = F()

class G(F):
    identity = A.element('e')
g = G()

def test_Monoid():
    # Monoid's operator must have identity
    raises(TypeError, lambda: Monoid('M', (A,), (f,)))

    M = Monoid('M', (A,), (g,))
    M_op = M.operator
    a, b, e = [A.element(i) for i in 'abe']

    assert M_op(a, e, evaluate=True) == M_op(e, a, evaluate=True) == a
    assert M_op(M_op(a, e), b, evaluate=True).args == (M_op, (a, b))
