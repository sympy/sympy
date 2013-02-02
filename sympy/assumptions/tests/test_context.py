from __future__ import with_statement
from sympy.assumptions import ask, Q
from sympy.assumptions.assume import assume
from sympy.abc import x, y

def test_assume():
    with assume(Q.integer(x)):
        assert ask(Q.integer(x))
    assert not ask(Q.integer(x))

def test_assume_nested():
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))
    with assume(Q.integer(x)):
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
        with assume(Q.integer(y)):
            assert ask(Q.integer(x))
            assert ask(Q.integer(y))
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))

def test_finally():
    try:
        with assume(Q.integer(x)):
            1/0
    except ZeroDivisionError:
        pass
    assert not ask(Q.integer(x))
