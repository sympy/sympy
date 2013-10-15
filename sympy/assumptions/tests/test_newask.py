from sympy.assumptions.newask import newask

from sympy import symbols, Q, assuming, Implies, MatrixSymbol

from sympy.utilities.pytest import raises, XFAIL

x, y = symbols('x y')

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

    assert newask(Q.zero(x), Q.nonzero(x)) is False
    assert newask(Q.positive(x), Q.zero(x)) is False
    assert newask(Q.real(x), Q.zero(x)) is True
    assert newask(Q.zero(x), Q.zero(x*y)) is None
    assert newask(Q.zero(x*y), Q.zero(x))

def test_zero():
    """
    Everything in this test doesn't work with ask, and most things would be
    very difficult or impossible to make work under the current handlers
    model.
    """
    assert newask(Q.zero(x) | Q.zero(y), Q.zero(x*y)) is True
    assert newask(Q.zero(x*y), Q.zero(x) | Q.zero(y)) is True

    assert newask(Implies(Q.zero(x), Q.zero(x*y))) is True

    # This one in particular requires computing the fixed-point of the
    # relevant facts, because going from Q.nonzero(x*y) -> ~Q.zero(x*y) and
    # Q.zero(x*y) -> Equivalent(Q.zero(x*y), Q.zero(x) | Q.zero(y)) takes two
    # steps.
    assert newask(Q.zero(x) | Q.zero(y), Q.nonzero(x*y)) is False

    assert newask(Q.zero(x), Q.zero(x**2)) is True

def test_zero_positive():
    assert newask(Q.zero(x + y), Q.positive(x) & Q.positive(y)) is False
    assert newask(Q.positive(x) & Q.positive(y), Q.zero(x + y)) is False
    assert newask(Q.nonzero(x + y), Q.positive(x) & Q.positive(y)) is True
    assert newask(Q.positive(x) & Q.positive(y), Q.nonzero(x + y)) is None

    # This one requires several levels of forward chaining
    assert newask(Q.zero(x*(x + y)), Q.positive(x) & Q.positive(y)) is False

def test_zero_pow():
    assert newask(Q.zero(x**y), Q.zero(x) & Q.positive(y)) is True
    assert newask(Q.zero(x**y), Q.nonzero(x) & Q.zero(y)) is False

    assert newask(Q.zero(x), Q.zero(x**y)) is True

    assert newask(Q.zero(x**y), Q.zero(x)) is None

@XFAIL
# Requires correct Q.square calculation first
def test_invertible():
    A = MatrixSymbol('A', 5, 5)
    B = MatrixSymbol('B', 5, 5)
    assert newask(Q.invertible(A*B), Q.invertible(A) & Q.invertible(B)) is True
    assert newask(Q.invertible(A), Q.invertible(A*B))
    assert newask(Q.invertible(A) & Q.invertible(B), Q.invertible(A*B))
