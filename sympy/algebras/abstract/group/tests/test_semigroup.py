from sympy import Set, Map, BinaryOperator
from sympy.algebras import Semigroup
from sympy.testing.pytest import raises

A = Set('A')

f = Map('f', domain=A**2, codomain=A)

class G(BinaryOperator):
    domain = A*A
    codomain = A
    associative = True
g = G()

def test_Semigroup():
    # Semigroup's operator must be associative
    raises(TypeError, lambda: Semigroup('SG', (A,), (f,)))

    SG = Semigroup('SG', (A,), (g,))
    SG_op = SG.operator
    a, b, c = [A.element(i) for i in 'abc']

    assert SG_op(a, SG_op(b, c), evaluate=True).args == (SG_op, (a, b, c))
    assert SG_op(SG_op(a, b), c, evaluate=True).args == (SG_op, (a, b, c))
