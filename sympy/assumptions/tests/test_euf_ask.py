# sympy/assumptions/tests/test_euf_propagation.py

from sympy.assumptions.euf_ask import euf_ask
from sympy.assumptions.ask import Q
from sympy import symbols

# simple symbols
x, y, z, u, v = symbols("x y z u v")

def test_direct_propagation():
    # Direct one-step propagation: same predicate on y
    assert euf_ask(Q.prime(x), Q.prime(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.even(x), Q.even(y) & Q.eq(x, y)) is True
    assert euf_ask(Q.positive(x), Q.positive(y) & Q.eq(x, y)) is True

def test_chain_propagation_two_steps():
    # x should inherit from y through z
    assert euf_ask(Q.prime(x), Q.prime(z) & Q.eq(z, y) & Q.eq(y, x)) is True
    assert euf_ask(Q.integer(x), Q.integer(z) & Q.eq(z, y) & Q.eq(y, x)) is True
    assert euf_ask(Q.real(x), Q.real(z) & Q.eq(z, y) & Q.eq(y, x)) is True

def test_chain_three_steps_various_preds():
    # Multiple-hop equality chain
    assumptions = Q.eq(x, y) & Q.eq(y, z) & Q.eq(z, u) & Q.eq(u, v)
    assert euf_ask(Q.rational(x), Q.rational(v) & assumptions) is True
    assert euf_ask(Q.algebraic(x), Q.algebraic(v) & assumptions) is True
    assert euf_ask(Q.transcendental(x), Q.transcendental(v) & assumptions) is True

def test_mixed_predicates_all_true():
    # x should be integer because v is integer through chain
    chain = Q.eq(x, y) & Q.eq(y, z) & Q.eq(z, v)
    assert euf_ask(Q.integer(x), Q.integer(v) & chain) is True
    # And positive because y is positive
    assert euf_ask(Q.positive(x), Q.positive(y) & chain) is True

def test_nested_conjunction_in_assumptions():
    # Nested & inside parentheses
    nested = (Q.prime(y) & Q.eq(y, z)) & (Q.eq(z, x))
    assert euf_ask(Q.prime(x), nested) is True
    chain2 = (Q.positive(y) & Q.eq(y, z)) & (Q.positive(z) & Q.eq(z, x))
    assert euf_ask(Q.positive(x), chain2) is True

def test_predicate_switch_symmetry():
    # Swapping x and y still propagates
    assert euf_ask(Q.prime(y), Q.prime(x) & Q.eq(x, y)) is True
    assert euf_ask(Q.odd(y), Q.odd(x) & Q.eq(x, y)) is True

def test_multiple_predicates_different_vars():
    # x inherits prime from z, inherits positive from y
    mixed = Q.prime(z) & Q.eq(z, x) & Q.positive(y) & Q.eq(y, x)
    assert euf_ask(Q.prime(x), mixed) is True
    assert euf_ask(Q.positive(x), mixed) is True
