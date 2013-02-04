from __future__ import with_statement
from sympy.assumptions import ask, Q
from sympy.assumptions.assume import assuming
from sympy.abc import x, y

def test_assuming():
    with assuming(Q.integer(x)):
        assert ask(Q.integer(x))
    assert not ask(Q.integer(x))

def test_assuming_nested():
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))
    with assuming(Q.integer(x)):
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
        with assuming(Q.integer(y)):
            assert ask(Q.integer(x))
            assert ask(Q.integer(y))
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))

def test_finally():
    try:
        with assuming(Q.integer(x)):
            1/0
    except ZeroDivisionError:
        pass
    assert not ask(Q.integer(x))
