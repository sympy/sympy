from sympy.assumptions.euf import euf_ask
from sympy.assumptions.ask import Q
from sympy.core import symbols

x, y, z= symbols("x y z")

def test_euf_basic_predicates():
    # Number theory: parity and integer
    assert euf_ask(Q.integer(x), Q.even(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.odd(x), Q.even(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.even(x), Q.odd(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.integer(x), Q.odd(y) & Q.eq(x, y)) is True

def test_euf_composite_and_prime():
    # Composite, prime, and propagation
    assert euf_ask(Q.prime(x), Q.composite(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.composite(x), Q.prime(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.integer(x), Q.composite(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.positive(x), Q.composite(y) & Q.eq(x, y)) is True

def test_euf_rational_irrational_real():
    # Rational, irrational, real
    assert euf_ask(Q.real(x), Q.rational(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.irrational(x), Q.rational(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.real(x), Q.irrational(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.rational(x), Q.irrational(y) & Q.eq(x, y)) is False

def test_euf_matrix_properties():
    # Matrix facts: orthogonal implies unitary, invertible, etc.
    assert euf_ask(Q.unitary(x), Q.orthogonal(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.positive_definite(x), Q.orthogonal(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.invertible(x), Q.orthogonal(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.fullrank(x), Q.orthogonal(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.square(x), Q.orthogonal(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.singular(x), Q.orthogonal(y) & Q.eq(x, y)) is False

def test_euf_exclusivity_and_contradiction():
    # Exclusivity: even/odd, contradiction detection
    assert euf_ask(Q.odd(x), Q.even(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.even(x), Q.odd(y) & Q.eq(x, y)) is False

def test_euf_equality_and_disequality():
    # Equality and disequality
    assert euf_ask(Q.eq(x, y), Q.eq(x, y)) is True
    assert euf_ask(Q.ne(x, y), Q.eq(x, y)) is False
    assert euf_ask(Q.eq(x, y), Q.ne(x, y)) is False
    assert euf_ask(Q.ne(x, y), Q.ne(x, y)) is True

def test_euf_ambiguous_cases():
    # Ambiguous/unknown: even + prime could be 2 (even+prime) or 5 (odd+prime)
    assert euf_ask(Q.even(x), Q.prime(y) & Q.eq(x, y)) is None
    # No info about x
    assert euf_ask(Q.positive(x), True) is None

def test_euf_matrix_exclusivity():
    # Matrix exclusivity: singular and invertible
    assert euf_ask(Q.singular(x), Q.invertible(x)) is False
    assert euf_ask(Q.invertible(x), Q.singular(x)) is False

def test_euf_reflexivity_and_transitivity():
    # Reflexivity and transitivity for equality
    assert euf_ask(Q.eq(x, x)) is True
    assert euf_ask(Q.eq(y, x), Q.eq(x, y)) is True
    # Transitivity: x=y, y=z => x=z
    assert euf_ask(Q.eq(x, z), Q.eq(x, y) & Q.eq(y, z)) is True
