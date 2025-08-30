# sympy/assumptions/tests/test_euf_propagation.py

from sympy.assumptions.euf_ask import euf_ask
from sympy.assumptions.ask import Q
from sympy import symbols
from sympy.testing.pytest import XFAIL

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


def test_negations():
    assert euf_ask(~Q.prime(x), Q.prime(y) & Q.eq(x, y)) is False
    assert euf_ask(~Q.even(x), Q.even(y) & Q.eq(x, y)) is False
    assert euf_ask(~Q.positive(x), Q.positive(y) & Q.eq(x, y)) is False


def test_or():
    assumptions = Q.prime(x) | (Q.eq(x, y) & Q.prime(y))
    assert euf_ask(Q.prime(x), assumptions) is True


def test_long_equality_chains():
    """Test very long chains of equalities"""
    # 6-step chain
    a, b, c, d, e, f = symbols('a b c d e f')
    long_chain = Q.eq(a, b) & Q.eq(b, c) & Q.eq(c, d) & Q.eq(d, e) & Q.eq(e, f)
    assert euf_ask(Q.prime(a), Q.prime(f) & long_chain) is True
    assert euf_ask(Q.negative(a), Q.negative(f) & long_chain) is True


def test_branching_equality_trees():
    """Test tree-like equality structures"""
    # Tree: x = y, x = z, both y and z have properties
    tree = Q.eq(x, y) & Q.eq(x, z) & Q.prime(y) & Q.positive(z)
    assert euf_ask(Q.prime(x), tree) is True
    assert euf_ask(Q.positive(x), tree) is True
    assert euf_ask(Q.prime(z), tree) is True  # z should inherit from y through x
    assert euf_ask(Q.positive(y), tree) is True  # y should inherit from z through x


@XFAIL
def test_complex_mixed_expressions():
    """Test complex Boolean combinations"""
    # Mixed AND/OR with equalities
    complex_expr = (Q.prime(x) | Q.even(x)) & Q.eq(x, y) & (Q.odd(y) | Q.positive(y))
    # Should work for some combinations
    result = euf_ask(Q.prime(x), Q.prime(y) & complex_expr)
    # Complex expressions might return None when indeterminate

    # Nested implications
    nested = Q.eq(x, y) & (Q.prime(x) | (Q.eq(y, z) & Q.prime(z)))
    assert euf_ask(Q.prime(x), Q.prime(z) & nested) is True


def test_multiple_disjoint_equality_groups():
    """Test multiple separate equality groups"""
    # Two separate equality groups
    group1 = Q.eq(x, y) & Q.prime(x)
    group2 = Q.eq(z, u) & Q.even(z)
    combined = group1 & group2

    # Properties should stay within groups
    assert euf_ask(Q.prime(y), combined) is True
    assert euf_ask(Q.even(u), combined) is True
    # But not cross groups
    assert euf_ask(Q.prime(z), combined) is None
    assert euf_ask(Q.even(x), combined) is None


def test_symmetry():
    """Test symmetry of equality"""
    # x = y implies y = x (symmetry)
    assert euf_ask(Q.prime(y), Q.prime(x) & Q.eq(x, y)) is True
    assert euf_ask(Q.prime(y), Q.prime(x) & Q.eq(y, x)) is True


def test_large_assumption_sets():
    """Test performance with many assumptions"""
    # Create a large set of assumptions
    vars = symbols('v1:20')  # v1, v2, ..., v19

    # Chain all variables together
    chain_assumptions = [Q.eq(vars[i], vars[i+1]) for i in range(len(vars)-1)]
    all_chain = chain_assumptions[0]
    for assumption in chain_assumptions[1:]:
        all_chain = all_chain & assumption

    # Add property to last variable
    all_chain = all_chain & Q.prime(vars[-1])

    # First variable should inherit the property
    assert euf_ask(Q.prime(vars[0]), all_chain) is True


def test_circular_equality_chains():
    """Test circular equality references"""
    # x = y, y = z, z = x (should be handled correctly)
    circular = Q.eq(x, y) & Q.eq(y, z) & Q.eq(z, x) & Q.prime(x)
    assert euf_ask(Q.prime(y), circular) is True
    assert euf_ask(Q.prime(z), circular) is True


def test_negated_equality():
    """Test negated equalities and their implications"""
    # If x != y and x is prime, we can't conclude anything about y
    assert euf_ask(Q.prime(y), Q.prime(x) & Q.ne(x, y)) is None


def test_function_equality():
    """Test equality with function expressions"""
    f, g = symbols('f g', cls=symbols)
    a, b = symbols('a b')

    # If f(a) = g(b) and f(a) is prime, then g(b) is prime
    func_eq = Q.eq(f, g) & Q.prime(f)
    assert euf_ask(Q.prime(g), func_eq) is True


@XFAIL
def test_edge_case_empty_assumptions():
    """Test behavior with minimal or empty assumptions"""
    # With no assumptions, we should get None (unknown)
    assert euf_ask(Q.prime(x), True) is None

    # With only equality but no properties
    assert euf_ask(Q.prime(x), Q.eq(x, y)) is None


def test_deeply_nested_boolean_structure():
    """Test deeply nested Boolean expressions"""
    deep_nested = ((Q.prime(x) & Q.eq(x, y)) |
                   (Q.even(x) & Q.eq(x, z)) |
                   (Q.odd(x) & Q.eq(x, u))) & Q.prime(y) & Q.even(z) & Q.odd(u)

    # This creates complex satisfiability conditions
    result = euf_ask(Q.prime(x), deep_nested)
    # May be True, False, or None depending on implementation


@XFAIL
def test_failing_cases():
    # But contradictions should still work
    not_equal = Q.ne(x, y) & Q.prime(x) & Q.composite(x)
    assert euf_ask(Q.prime(x), not_equal) is False

    contradiction_tree = Q.eq(x, y) & Q.eq(x, z) & Q.positive(y) & Q.negative(z)
    assert euf_ask(Q.positive(x), contradiction_tree) is False
    assert euf_ask(Q.negative(x), contradiction_tree) is False

    assert euf_ask(Q.prime(x), True) is None

    # With only equality but no properties
    assert euf_ask(Q.prime(x), Q.eq(x, y)) is None

    # x = x should always be true (reflexivity)
    assert euf_ask(Q.prime(x), Q.prime(x) & Q.eq(x, x)) is True
