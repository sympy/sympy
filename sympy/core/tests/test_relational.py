from sympy import symbols, oo
from sympy.core.relational import Relational, Equality, StrictInequality, \
    Rel, Eq, Lt, Le, Gt, Ge, Ne

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


def test_wrappers():
    e = x+x**2

    res = Relational(y, e, '==')
    assert Rel(y, x+x**2, '==') == res
    assert Eq(y, x+x**2) == res

    res = Relational(y, e, '<')
    assert Lt(y, x+x**2) == res

    res = Relational(y, e, '<=')
    assert Le(y, x+x**2) == res

    res = Relational(y, e, '>')
    assert Gt(y, x+x**2) == res

    res = Relational(y, e, '>=')
    assert Ge(y, x+x**2) == res

    res = Relational(y, e, '!=')
    assert Ne(y, x+x**2) == res

def test_Eq():

    assert Eq(x**2) == Eq(x**2, 0)
    assert Eq(x**2) != Eq(x**2, 1)

def test_rel_Infinity():
    assert (oo > oo) is False
    assert (oo > -oo) is True
    assert (oo > 1) is True
    assert (oo < oo) is False
    assert (oo < -oo) is False
    assert (oo < 1) is False
    assert (oo >= oo) is True
    assert (oo >= -oo) is True
    assert (oo >= 1) is True
    assert (oo <= oo) is True
    assert (oo <= -oo) is False
    assert (oo <= 1) is False
    assert (-oo > oo) is False
    assert (-oo > -oo) is False
    assert (-oo > 1) is False
    assert (-oo < oo) is True
    assert (-oo < -oo) is False
    assert (-oo < 1) is True
    assert (-oo >= oo) is False
    assert (-oo >= -oo) is True
    assert (-oo >= 1) is False
    assert (-oo <= oo) is True
    assert (-oo <= -oo) is True
    assert (-oo <= 1) is True
