from __future__ import annotations
from sympy.assumptions.satask import satask
from sympy.assumptions.ask import Q, ask

from sympy.core import symbols, Symbol
from sympy.functions.elementary.complexes import Abs
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.core.numbers import I, oo

from sympy.testing.pytest import raises, XFAIL
x, y, z = symbols("x y z", real=True)

def test_satask_lra():
    im = Symbol('im', imaginary=True)

    # test preprocessing of unequalities is working correctly
    assert satask(Q.eq(x, 1), ~Q.ne(x, 0)) is False
    assert satask(Q.eq(x, 0), ~Q.ne(x, 0)) is True
    assert satask(~Q.ne(x, 0), Q.eq(x, 0)) is True
    assert satask(~Q.eq(x, 0), Q.eq(x, 0)) is False
    assert satask(Q.ne(x, 0), Q.eq(x, 0)) is False

    # basic tests
    assert satask(Q.ne(x, x)) is False
    assert satask(Q.eq(x, x)) is True
    assert satask(Q.gt(x, 0), Q.gt(x, 1)) is True

    # check that True/False are handled
    assert satask(Q.gt(x, 0), True) is None
    raises(ValueError, lambda: satask(Q.gt(x, 0), False))

    # check imaginary numbers are correctly handled
    # (im * I).is_real returns True so this is an edge case
    assert satask(Q.gt(im * I, 0), Q.gt(im * I, 1)) is None

    # check matrix inputs
    X = MatrixSymbol("X", 2, 2)
    assert satask(Q.lt(X, 2) & Q.gt(X, 3)) is None


def test_old_assumptions():
    # test unhandled old assumptions
    w = symbols("w")
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", rational=False, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", odd=True, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", even=True, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", prime=True, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", composite=True, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", integer=True, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None
    w = symbols("w", integer=False, real=True)
    assert satask(Q.lt(w, 2) & Q.gt(w, 3)) is None

    # test handled
    w = symbols("w", positive=True, real=True)
    assert satask(Q.le(w, 0)) is False
    assert satask(Q.gt(w, 0)) is True
    w = symbols("w", negative=True, real=True)
    assert satask(Q.lt(w, 0)) is True
    assert satask(Q.ge(w, 0)) is False
    w = symbols("w", zero=True, real=True)
    assert satask(Q.eq(w, 0)) is True
    assert satask(Q.ne(w, 0)) is False
    w = symbols("w", nonzero=True, real=True)
    assert satask(Q.ne(w, 0)) is True
    assert satask(Q.eq(w, 1)) is None
    w = symbols("w", nonpositive=True, real=True)
    assert satask(Q.le(w, 0)) is True
    assert satask(Q.gt(w, 0)) is False
    w = symbols("w", nonnegative=True, real=True)
    assert satask(Q.ge(w, 0)) is True
    assert satask(Q.lt(w, 0)) is False


def test_new_assumptions():
    # Every expression has to be known to be real for the theory solver to
    # be able to answer, and here that is what the new assumptions say.
    a, b, c = symbols("a b c")
    real = Q.real(a) & Q.real(b) & Q.real(c)

    # Q.real itself
    assert satask(a > c, (a > b) & (b > c) & real) is True
    assert satask(a <= c, (a <= b) & (b <= c) & real) is True
    assert satask(Q.eq(a, c), Q.eq(a, b) & Q.eq(b, c) & real) is True
    assert satask(Q.ne(a, c), (a < b) & (b < c) & real) is True
    assert satask(a + b > 2, (a < 0) & (b < 0) & real) is False
    assert satask(Q.lt(a, 2) & Q.gt(a, 3), Q.real(a)) is False
    assert satask(a > 0, (Q.gt(a, 1) | Q.gt(a, 2)) & Q.real(a)) is True

    # a predicate that implies Q.real is just as good
    assert satask(a - b > 0, Q.nonnegative(a) & Q.negative(b)) is True
    assert satask(Q.ne(a, b), Q.nonnegative(a) & Q.negative(b)) is True
    assert satask(a > b, Q.positive(a) & Q.zero(b)) is True
    assert satask(a > 0, Q.nonnegative(a) & Q.nonzero(a)) is True
    assert satask(a < 0, Q.negative(a)) is True
    assert satask(Q.extended_positive(a), Q.gt(a, 2) & Q.real(a)) is True
    assert satask(Q.positive_infinite(a), Q.real(a)) is False

    # the old and the new assumptions can be mixed
    assert satask(x > c, (x > b) & (b > c) & Q.positive(b) & Q.real(c)) is True

    # an expression is real when the expressions it is built from are
    assert satask(2*a > a, Q.positive(a)) is True
    assert satask(a/2 > 0, Q.positive(a)) is True
    assert satask(a + 1 > a, Q.real(a)) is True
    assert satask(Abs(a) >= 0, Q.real(a)) is True


def test_unhandled_new_assumptions():
    a, b, c = symbols("a b c")
    real = Q.real(a) & Q.real(b) & Q.real(c)

    # An inequality only says that its expressions are extended real, and
    # extended real arithmetic is not what the theory solver implements.
    assert satask(a > c, (a > b) & (b > c)) is None
    assert satask(a + 1 > a) is None

    # For these, None is not just the answer that is given but the answer
    # that is right: a could be oo, and both 2*oo > oo and oo + 1 > oo are
    # False, so nothing stronger may be concluded.
    assert satask(2*a > a, a > 0) is None
    assert satask(2*a > a, Q.extended_positive(a)) is None
    assert satask(a + 1 > a, Q.extended_real(a)) is None

    # what makes the expressions real has to be assumed, not proposed
    assert satask(Q.lt(a, 2) & Q.gt(a, 3) & Q.real(a)) is None

    # predicates the theory solver cannot account for
    assert satask(a > 0, Q.gt(a, 1) & Q.integer(a)) is None
    assert satask(Q.gt(a, 0) & Q.lt(a, 1), Q.prime(a)) is None
    assert satask(a > 0, Q.gt(a, 1) & Q.rational(a)) is None
    assert satask(a > 0, Q.gt(a, 1) & Q.irrational(a)) is None

    # a predicate that leaves the expressions complex
    assert satask(a > 0, Q.gt(a, 1) & Q.complex(a)) is None

    # nonlinear arithmetic
    assert satask(a*c <= b*c, (a <= b) & (c > 0) & real) is None

    # infinity
    assert satask(Q.lt(a, oo), Q.real(a)) is None


def test_rel_queries():
    assert ask(Q.lt(x, 2) & Q.gt(x, 3)) is False
    assert ask(Q.positive(x - z), (x > y) & (y > z)) is True
    assert ask(x + y > 2, (x < 0) & (y <0)) is False
    assert ask(x > z, (x > y) & (y > z)) is True


def test_unhandled_queries():
    X = MatrixSymbol("X", 2, 2)
    assert ask(Q.lt(X, 2) & Q.gt(X, 3)) is None


def test_all_pred():
    # test usable pred
    assert satask(Q.extended_positive(x), (x > 2)) is True
    assert satask(Q.positive_infinite(x)) is False
    assert satask(Q.negative_infinite(x)) is False

    # test disallowed pred
    assert satask((x > 0), (x > 2) & Q.prime(x)) is None
    assert satask((x > 0), (x > 2) & Q.composite(x)) is None
    assert satask((x > 0), (x > 2) & Q.odd(x)) is None
    assert satask((x > 0), (x > 2) & Q.even(x)) is None
    assert satask((x > 0), (x > 2) & Q.integer(x)) is None


def test_number_line_properties():
    # From:
    # https://en.wikipedia.org/wiki/Inequality_(mathematics)#Properties_on_the_number_line

    a, b, c = symbols("a b c", real=True)

    # Transitivity
    # If a <= b and b <= c, then a <= c.
    assert ask(a <= c, (a <= b) & (b <= c)) is True
    # If a <= b and b < c, then a < c.
    assert ask(a < c, (a <= b) & (b < c)) is True
    # If a < b and b <= c, then a < c.
    assert ask(a < c, (a < b) & (b <= c)) is True

    # Addition and subtraction
    # If a <= b, then a + c <= b + c and a - c <= b - c.
    assert ask(a + c <= b + c, a <= b) is True
    assert ask(a - c <= b - c, a <= b) is True


@XFAIL
def test_failing_number_line_properties():
    # From:
    # https://en.wikipedia.org/wiki/Inequality_(mathematics)#Properties_on_the_number_line

    a, b, c = symbols("a b c", real=True)

    # Multiplication and division
    # If a <= b and c > 0, then ac <= bc and a/c <= b/c. (True for non-zero c)
    assert ask(a*c <= b*c, (a <= b) & (c > 0) & ~ Q.zero(c)) is True
    assert ask(a/c <= b/c, (a <= b) & (c > 0) & ~ Q.zero(c)) is True
    # If a <= b and c < 0, then ac >= bc and a/c >= b/c. (True for non-zero c)
    assert ask(a*c >= b*c, (a <= b) & (c < 0) & ~ Q.zero(c)) is True
    assert ask(a/c >= b/c, (a <= b) & (c < 0) & ~ Q.zero(c)) is True

    # Additive inverse
    # If a <= b, then -a >= -b.
    assert ask(-a >= -b, a <= b) is True

    # Multiplicative inverse
    # For a, b that are both negative or both positive:
    # If a <= b, then 1/a >= 1/b .
    assert ask(1/a >= 1/b, (a <= b) & Q.positive(x) & Q.positive(b)) is True
    assert ask(1/a >= 1/b, (a <= b) & Q.negative(x) & Q.negative(b)) is True


def test_equality():
    # test symmetry and reflexivity
    assert ask(Q.eq(x, x)) is True
    assert ask(Q.eq(y, x), Q.eq(x, y)) is True
    assert ask(Q.eq(y, x), ~Q.eq(z, z) | Q.eq(x, y)) is True

    # test transitivity
    assert ask(Q.eq(x,z), Q.eq(x,y) & Q.eq(y,z)) is True


@XFAIL
def test_equality_failing():
    # Note that implementing the substitution property of equality
    # most likely requires a redesign of the new assumptions.
    # See issue #25485 for why this is the case and general ideas
    # about how things could be redesigned.

    # test substitution property
    assert ask(Q.prime(x), Q.eq(x, y) & Q.prime(y)) is True
    assert ask(Q.real(x), Q.eq(x, y) & Q.real(y)) is True
    assert ask(Q.imaginary(x), Q.eq(x, y) & Q.imaginary(y)) is True
