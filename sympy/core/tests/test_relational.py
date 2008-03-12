from sympy import symbols
from sympy.core.relational import Relational, Equality, StrictInequality

x,y,z = symbols('xyz')


def test_rel_ne():
    Relational(x, y, '!=')  # this used to raise


def test_rel_subs():
    e = Relational(x, y, '==')
    e = e.subs(x,z)

    assert isinstance(e, Equality)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '<')
    e = e.subs(x,z)

    assert isinstance(e, StrictInequality)
    assert e.lhs == z
    assert e.rhs == y

