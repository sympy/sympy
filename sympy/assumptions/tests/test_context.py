from sympy.assumptions import ask, Q
from sympy.assumptions.assume import assume
from sympy.abc import x

def test_context_manager():
    with assume(Q.integer(x)):
        assert ask(Q.integer(x))
    assert not ask(Q.integer(x))
