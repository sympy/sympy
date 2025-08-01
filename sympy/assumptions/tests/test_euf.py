from sympy.assumptions.euf import euf_ask
from sympy.assumptions.ask import Q
from sympy.core import symbols
from sympy.matrices.expressions.matexpr import MatrixSymbol

x, y, z= symbols("x y z")
M = MatrixSymbol("M", 3, 3)

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


# Algebraic and Transcendental Numbers
def test_algebraic_transcendental():
    assert euf_ask(Q.algebraic(x), Q.algebraic(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.transcendental(x), Q.algebraic(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.complex(x), Q.transcendental(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.integer(x), Q.transcendental(y) & Q.eq(x, y)) is False


# Matrix Structural Properties
def test_matrix_structural_properties():
    assert euf_ask(Q.diagonal(M), Q.diagonal(M)) is True
    assert euf_ask(Q.upper_triangular(M), Q.diagonal(M)) is True
    assert euf_ask(Q.lower_triangular(M), Q.diagonal(M)) is True
    assert euf_ask(Q.symmetric(M), Q.diagonal(M)) is True
    assert euf_ask(Q.normal(M), Q.unitary(M)) is True  # Unitary implies normal


# Extended Real Number System
def test_extended_real_properties():
    assert euf_ask(Q.extended_positive(x), Q.positive_infinite(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.finite(x), Q.positive_infinite(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.extended_nonzero(x), Q.negative_infinite(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.real(x), Q.positive_infinite(y) & Q.eq(x, y)) is False


# Zero and Special Number Properties
def test_zero_properties():
    assert euf_ask(Q.zero(x), Q.zero(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.nonzero(x), Q.zero(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.even(x), Q.zero(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.prime(x), Q.zero(y) & Q.eq(x, y)) is False


# Triangular Matrix Properties
def test_triangular_matrices():
    assert euf_ask(Q.triangular(M), Q.lower_triangular(M)) is True
    assert euf_ask(Q.upper_triangular(M), Q.unit_triangular(M)) is None


# Complex Number Components
def test_complex_components():
    assert euf_ask(Q.real_elements(M), Q.integer_elements(M)) is True
    assert euf_ask(Q.complex_elements(M), Q.real_elements(M)) is True


# Matrix Rank Properties
def test_matrix_rank_properties():
    assert euf_ask(Q.singular(M), Q.orthogonal(M)) is False
    assert euf_ask(Q.invertible(M), Q.unitary(M)) is True


# Number Type Relationships
def test_number_type_relationships():
    assert euf_ask(Q.irrational(x), Q.rational(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.real(x), Q.irrational(y) & Q.eq(x, y)) is True


# Prime Number Exclusive Properties
def test_prime_exclusivity():
    assert euf_ask(Q.composite(x), Q.prime(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.even(x), Q.prime(y) & Q.eq(x, y)) is None  # Could be 2
    assert euf_ask(Q.odd(x), Q.prime(y) & Q.eq(x, y)) is None  # Except 2


# Extended Real Comparisons
def test_extended_real_inequalities():
    assert euf_ask(Q.extended_nonnegative(x), Q.zero(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.extended_positive(x), Q.positive_infinite(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.extended_nonpositive(x), Q.negative_infinite(y) & Q.eq(x, y)) is True


# Matrix Operation Properties
def test_matrix_operations():
    assert euf_ask(Q.positive_definite(M), Q.orthogonal(M)) is True
    assert euf_ask(Q.normal(M), Q.unitary(M)) is True
    assert euf_ask(Q.square(M), Q.diagonal(M)) is True


# Special Number Sets
def test_special_number_sets():
    assert euf_ask(Q.finite(x), Q.algebraic(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.infinite(x), Q.transcendental(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.commutative(x), Q.prime(y) & Q.eq(x, y)) is True


# Complex Real Interactions
def test_complex_real_interactions():
    assert euf_ask(Q.complex(x), Q.imaginary(y) & Q.eq(x, y)) is True


# Integer Type Properties
def test_integer_type_properties():
    assert euf_ask(Q.nonnegative(x), Q.prime(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.positive(x), Q.composite(y) & Q.eq(x, y)) is True


# Matrix Special Properties
def test_matrix_special_properties():
    assert euf_ask(Q.triangular(M), Q.unit_triangular(M)) is True
    assert euf_ask(Q.symmetric(M), Q.normal(M)) is None  # Not all normal matrices are symmetric


# Number Line Extremes
def test_number_line_extremes():
    assert euf_ask(Q.finite(x), Q.negative_infinite(y) & Q.eq(x, y)) is False
    assert euf_ask(Q.extended_nonzero(x), Q.positive_infinite(y) & Q.eq(x, y)) is True


# Composite Number Implications
def test_composite_implications():
    assert euf_ask(Q.positive(x), Q.composite(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.transcendental(x), Q.composite(y) & Q.eq(x, y)) is False
