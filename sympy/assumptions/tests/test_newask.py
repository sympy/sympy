from sympy.assumptions.newask import newask

from sympy import symbols, Q, assuming, Implies, MatrixSymbol, I, pi, Rational

from sympy.utilities.pytest import raises, XFAIL

x, y, z = symbols('x y z')

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

    assert newask(Q.positive(pi*x*y + 1), Q.positive(x) & Q.positive(y)) is True
    assert newask(Q.positive(pi*x*y - 5), Q.positive(x) & Q.positive(y)) is None

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

def test_prime():
    assert newask(Q.prime(5)) is True
    assert newask(Q.prime(6)) is False
    assert newask(Q.prime(-5)) is False

    assert newask(Q.prime(x*y), Q.integer(x) & Q.integer(y)) is None
    assert newask(Q.prime(x*y), Q.prime(x) & Q.prime(y)) is False

def test_old_assump():
    assert newask(Q.positive(1)) is True
    assert newask(Q.positive(-1)) is False
    assert newask(Q.positive(0)) is False
    assert newask(Q.positive(I)) is False
    assert newask(Q.positive(pi)) is True

    assert newask(Q.negative(1)) is False
    assert newask(Q.negative(-1)) is True
    assert newask(Q.negative(0)) is False
    assert newask(Q.negative(I)) is False
    assert newask(Q.negative(pi)) is False

    assert newask(Q.zero(1)) is False
    assert newask(Q.zero(-1)) is False
    assert newask(Q.zero(0)) is True
    assert newask(Q.zero(I)) is False
    assert newask(Q.zero(pi)) is False

    assert newask(Q.nonzero(1)) is True
    assert newask(Q.nonzero(-1)) is True
    assert newask(Q.nonzero(0)) is False
    assert newask(Q.nonzero(I)) is False
    assert newask(Q.nonzero(pi)) is True

    assert newask(Q.nonpositive(1)) is False
    assert newask(Q.nonpositive(-1)) is True
    assert newask(Q.nonpositive(0)) is True
    assert newask(Q.nonpositive(I)) is False
    assert newask(Q.nonpositive(pi)) is False

    assert newask(Q.nonnegative(1)) is True
    assert newask(Q.nonnegative(-1)) is False
    assert newask(Q.nonnegative(0)) is True
    assert newask(Q.nonnegative(I)) is False
    assert newask(Q.nonnegative(pi)) is True

def test_irrational():
    assert newask(Q.irrational(2)) is False
    assert newask(Q.rational(2)) is True
    assert newask(Q.irrational(pi)) is True
    assert newask(Q.rational(pi)) is False
    assert newask(Q.irrational(I)) is False
    assert newask(Q.rational(I)) is False

    assert newask(Q.irrational(x*y*z), Q.irrational(x) & Q.irrational(y) &
        Q.rational(z)) is None
    assert newask(Q.irrational(x*y*z), Q.irrational(x) & Q.rational(y) &
        Q.rational(z)) is True
    assert newask(Q.irrational(pi*x*y), Q.rational(x) & Q.rational(y)) is True

    assert newask(Q.irrational(x + y + z), Q.irrational(x) & Q.irrational(y) &
        Q.rational(z)) is None
    assert newask(Q.irrational(x + y + z), Q.irrational(x) & Q.rational(y) &
        Q.rational(z)) is True
    assert newask(Q.irrational(pi + x + y), Q.rational(x) & Q.rational(y)) is True

    assert newask(Q.irrational(x*y*z), Q.rational(x) & Q.rational(y) &
        Q.rational(z)) is False
    assert newask(Q.rational(x*y*z), Q.rational(x) & Q.rational(y) &
        Q.rational(z)) is True

    assert newask(Q.irrational(x + y + z), Q.rational(x) & Q.rational(y) &
        Q.rational(z)) is False
    assert newask(Q.rational(x + y + z), Q.rational(x) & Q.rational(y) &
        Q.rational(z)) is True

def test_even():
    assert newask(Q.even(2)) is True
    assert newask(Q.even(3)) is False

    assert newask(Q.even(x*y), Q.even(x) & Q.odd(y)) is True
    assert newask(Q.even(x*y), Q.even(x) & Q.integer(y)) is True
    assert newask(Q.even(x*y), Q.even(x) & Q.even(y)) is True
    assert newask(Q.even(x*y), Q.odd(x) & Q.odd(y)) is False
    assert newask(Q.even(x*y), Q.even(x)) is None
    assert newask(Q.even(x*y), Q.odd(x) & Q.integer(y)) is None
    assert newask(Q.even(x*y), Q.odd(x) & Q.odd(y)) is False

    assert newask(Q.even(abs(x)), Q.even(x)) is True
    assert newask(Q.even(abs(x)), Q.odd(x)) is False
    assert newask(Q.even(x), Q.even(abs(x))) is None # x could be complex

def test_odd():
    assert newask(Q.odd(2)) is False
    assert newask(Q.odd(3)) is True

    assert newask(Q.odd(x*y), Q.even(x) & Q.odd(y)) is False
    assert newask(Q.odd(x*y), Q.even(x) & Q.integer(y)) is False
    assert newask(Q.odd(x*y), Q.even(x) & Q.even(y)) is False
    assert newask(Q.odd(x*y), Q.odd(x) & Q.odd(y)) is True
    assert newask(Q.odd(x*y), Q.even(x)) is None
    assert newask(Q.odd(x*y), Q.odd(x) & Q.integer(y)) is None
    assert newask(Q.odd(x*y), Q.odd(x) & Q.odd(y)) is True

    assert newask(Q.odd(abs(x)), Q.even(x)) is False
    assert newask(Q.odd(abs(x)), Q.odd(x)) is True
    assert newask(Q.odd(x), Q.odd(abs(x))) is None # x could be complex


def test_integer():
    assert newask(Q.integer(1)) is True
    assert newask(Q.integer(Rational(1, 2))) is False

    assert newask(Q.integer(x + y), Q.integer(x) & Q.integer(y)) is True
    assert newask(Q.integer(x + y), Q.integer(x)) is None

    assert newask(Q.integer(x + y), Q.integer(x) & ~Q.integer(y)) is False
    assert newask(Q.integer(x + y + z), Q.integer(x) & Q.integer(y) &
        ~Q.integer(z)) is False
    assert newask(Q.integer(x + y + z), Q.integer(x) & ~Q.integer(y) &
        ~Q.integer(z)) is None
    assert newask(Q.integer(x + y + z), Q.integer(x) & ~Q.integer(y)) is None
    assert newask(Q.integer(x + y), Q.integer(x) & Q.irrational(y)) is False

    assert newask(Q.integer(x*y), Q.integer(x) & Q.integer(y)) is True
    assert newask(Q.integer(x*y), Q.integer(x)) is None

    assert newask(Q.integer(x*y), Q.integer(x) & ~Q.integer(y)) is None
    assert newask(Q.integer(x*y), Q.integer(x) & ~Q.rational(y)) is False
    assert newask(Q.integer(x*y*z), Q.integer(x) & Q.integer(y) &
        ~Q.rational(z)) is False
    assert newask(Q.integer(x*y*z), Q.integer(x) & ~Q.rational(y) &
        ~Q.rational(z)) is None
    assert newask(Q.integer(x*y*z), Q.integer(x) & ~Q.rational(y)) is None
    assert newask(Q.integer(x*y), Q.integer(x) & Q.irrational(y)) is False

def test_abs():
    assert newask(Q.nonnegative(abs(x))) is True
    assert newask(Q.positive(abs(x)), ~Q.zero(x)) is True
    assert newask(Q.zero(x), ~Q.zero(abs(x))) is False
    assert newask(Q.zero(x), Q.zero(abs(x))) is True
    assert newask(Q.nonzero(x), ~Q.zero(abs(x))) is None # x could be complex
    assert newask(Q.zero(abs(x)), Q.zero(x)) is True
