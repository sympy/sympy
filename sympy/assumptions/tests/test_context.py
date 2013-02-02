from sympy.assumptions import ask, Q
from sympy.assumptions.assume import assume
from sympy.abc import x, y

def test_assume():
    with assume(Q.integer(x)):
        assert ask(Q.integer(x))
    assert not ask(Q.integer(x))

def test_assume_nested():
    with assume(Q.integer(x)):
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
        assert not ask(Q.integer(x + y))
        with assume(Q.integer(y)):
            assert ask(Q.integer(x))
            assert ask(Q.integer(y))
            assert ask(Q.integer(x + y))
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
        assert not ask(Q.integer(x + y))
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))
