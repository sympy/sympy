from sympy.assumptions import Q
from sympy.core.relational import Relational, Eq, Ne

from sympy.abc import x, y


def test_relational_instancecheck():
    assert isinstance(Q.eq(x,y), Relational)
    assert isinstance(Q.eq(x,y), Eq)
    assert isinstance(Q.ne(x,y), Relational)
    assert isinstance(Q.ne(x,y), Ne)
