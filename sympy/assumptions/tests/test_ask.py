from sympy import ask, Q, symbols, pi
import pytest


def test_ask_basic_numeric():
    # rational check
    assert ask(Q.rational(pi)) is False


def test_ask_even_integer():
    x, y = symbols('x y')

    # x even and y integer ⇒ x*y even
    assert ask(Q.even(x*y), Q.even(x) & Q.integer(y)) is True


def test_ask_prime_expression():
    x = symbols('x')

    # 4*x cannot be prime if x integer
    assert ask(Q.prime(4*x), Q.integer(x)) is False


def test_ask_unknown_case():
    x = symbols('x')

    # cannot determine parity without assumption
    assert ask(Q.odd(3*x)) is None


def test_ask_inconsistent_assumptions():
    x = symbols('x')

    with pytest.raises(ValueError):
        ask(Q.integer(x), Q.even(x) & Q.odd(x))


def test_ask_boolean_proposition_required():
    x = symbols('x')

    with pytest.raises(TypeError):
        ask(x)


def test_ask_boolean_assumptions_required():
    x = symbols('x')

    with pytest.raises(TypeError):
        ask(Q.integer(x), x)


def test_positive_assumption():
    x = symbols('x')

    assert ask(Q.positive(x), Q.positive(x)) is True


def test_negative_assumption():
    x = symbols('x')

    assert ask(Q.negative(x), Q.negative(x)) is True


def test_nonzero_from_positive():
    x = symbols('x')

    # positive ⇒ nonzero
    assert ask(Q.nonzero(x), Q.positive(x)) is True