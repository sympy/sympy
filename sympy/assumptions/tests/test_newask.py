from sympy.assumptions.newask import newask

from sympy import symbols, Q, assuming

from sympy.utilities.pytest import raises

x = symbols('x')

def test_newask():
    # No relevant facts
    assert newask(Q.real(x), Q.real(x)) is True
    assert newask(Q.real(x), ~Q.real(x)) is False
    assert newask(Q.real(x)) is None

    assert newask(Q.real(x), Q.positive(x)) is True
    assert newask(Q.positive(x), Q.real(x)) is None
    assert newask(Q.real(x), ~Q.positive(x)) is None
    assert newask(Q.positive(x), ~Q.real(x)) is False

    raises(ValueError, lambda: newask(Q.real(x), Q.real(x) & ~Q.real(x)))

    with assuming(Q.positive(x)):
        assert newask(Q.real(x)) is True
        assert newask(~Q.positive(x)) is False
        raises(ValueError, lambda: newask(Q.real(x), ~Q.positive(x)))
