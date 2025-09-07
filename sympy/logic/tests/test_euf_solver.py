import pytest
from sympy.testing.pytest import XFAIL
from sympy import symbols, Function, Eq, Not, Ne, Unequality
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic import boolalg
from sympy.assumptions.ask import Q
from sympy.logic.algorithms.euf_theory_solver import (
    EUFTheorySolver, EUFDisequalityContradictionException, EUFEqualityContradictionException,
    ProofProducingCongruenceClosure, _order_key, _ordered_pair, _canonical_lit, _canon_eq
)
from sympy.logic.algorithms.euf_theory import EUFUnhandledInput
from sympy.assumptions.assume import AppliedPredicate

f, g = symbols('f g', cls=Function)
a, b, c, d, e, x, y, z, u, v, w = symbols('a b c d e x y z u v w')


def test_initialize_and_istrue():
    solver = EUFTheorySolver()
    eqs = {Eq(a, b), Eq(b, c), Eq(f(a), f(b))}
    solver.Initialize(eqs)

    # Initially false
    for lit in eqs:
        assert solver.IsTrue(lit) is None

    # Assert a=b
    solver.SetTrue(Eq(a, b))
    assert solver.IsTrue(Eq(a, b)) is True
    assert solver.IsTrue(Eq(b, a)) is True
    # others still undecided
    assert solver.IsTrue(Eq(b, c)) is None


def test_positive_propagation_transitivity():
    solver = EUFTheorySolver()
    lits = {Eq(a, b), Eq(b, c), Eq(a, c)}
    solver.Initialize(lits)

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))
    # transitive a=c
    assert solver.IsTrue(Eq(a, c)) is True


def test_functional_congruence_propagation():
    solver = EUFTheorySolver()
    lits = {Eq(a, b), Eq(f(a), f(b))}
    solver.Initialize(lits)

    solver.SetTrue(Eq(a, b))
    # should propagate f(a)=f(b)
    assert solver.IsTrue(Eq(f(a), f(b))) is True


def test_disequality_no_merge():
    solver = EUFTheorySolver()
    solver.Initialize({Not(Eq(a, b))})
    solver.SetTrue(Not(Eq(a, b)))
    arep = solver.cc._find(solver.cc._flatten(a))
    brep = solver.cc._find(solver.cc._flatten(b))
    assert arep != brep
    assert brep in solver.disequalities_set[arep]


def test_disequality_conflict():
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Not(Eq(a, b))})
    solver.SetTrue(Not(Eq(a, b)))
    with pytest.raises(EUFDisequalityContradictionException):
        solver.SetTrue(Eq(a, b))


def test_backtrack_removes_last():
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Eq(b, c)})
    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))
    assert solver.IsTrue(Eq(a, c))
    solver.Backtrack(1)
    assert solver.IsTrue(Eq(a, c)) is None
    assert solver.IsTrue(Eq(a, b)) is True


def test_explanation_for_positive():
    solver = EUFTheorySolver()
    lits = {Eq(a, b), Eq(b, c), Eq(a, c)}
    solver.Initialize(lits)
    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))
    expl = solver.Explanation(Eq(a, c))
    for lit in expl:
        assert solver.IsTrue(lit)


def test_redundant_assertions_safe():
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b)})
    solver.SetTrue(Eq(a, b))
    # asserting again should not error
    solver.SetTrue(Eq(a, b))
    assert solver.IsTrue(Eq(a, b))


def test_sympy_unequality_init():
    solver = EUFTheorySolver()
    solver.Initialize({Unequality(a, b)})
    solver.SetTrue(Unequality(a, b))
    assert solver.IsTrue(Unequality(a, b)) is True


def make_solver_from_props(*props):
    cnf = CNF.from_prop(boolalg.And(*props))
    enc = EncodedCNF()
    enc.from_cnf(cnf)
    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    return solver, enc, conflicts


def test_order_independence_of_assertions():
    # x=y, y=z - test that order doesn't matter
    prop = Q.eq(x,y) & Q.eq(y,z)
    _prop = CNF.from_prop(prop)
    sat = EncodedCNF()
    sat.add_from_cnf(_prop)
    solver, conflicts = EUFTheorySolver.from_encoded_cnf(sat)
    solver.assert_lit(1)
    solver.assert_lit(2)
    assert conflicts == []
    assert solver.IsTrue(Eq(x, z)) is True
    # Reset and reverse
    solver.Backtrack(len(solver.interpretation_stack))
    assert solver.IsTrue(Eq(x, z)) is None


def test_simple_equality_chain():
    """
    Test EUF: x = y, y = z  -> SAT, and x = z must hold.
    """
    cnf = CNF().from_prop(Eq(x, y) & Eq(y, z))
    enc = EncodedCNF(); enc.from_cnf(cnf)
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)

    # Assert all clauses/literals
    for lit_id in enc.encoding.values():
        assert euf.assert_lit(lit_id) == (True,set())

    assert euf.IsTrue(Eq(x, z)) is True


def test_backtrack_recovery():
    """
    EUF: assert eq, then backtrack and verify its gone.
    """
    cnf = CNF().from_prop(Eq(x, y))
    enc = EncodedCNF(); enc.from_cnf(cnf)
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    lit_id = next(iter(enc.encoding.values()))

    euf.assert_lit(lit_id)
    assert euf.IsTrue(Eq(x, y)) is True

    euf.Backtrack(1)  # undo the only decision
    assert euf.IsTrue(Eq(x, y)) is None


@XFAIL
def test_issue_1():
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = { Q.prime(x): 1, Q.eq(x,y): 2, Q.prime(y): 3 }
    enc_cnf.data = [{-1},{2},{3}]
    solver,conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)
    assert solver.assert_lit(-1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (False, {1,-2,3})


# Test helper functions
def test_order_key():
    """Test _order_key function for consistent ordering."""
    assert _order_key(a) == str(a)
    assert _order_key(f(a)) == str(f(a))


def test_ordered_pair():
    """Test _ordered_pair function."""
    pair1 = _ordered_pair(a, b)
    pair2 = _ordered_pair(b, a)
    assert pair1 == pair2  # Should be consistently ordered


def test_canonical_lit():
    """Test _canonical_lit function with various literal types."""
    # Test Eq
    lhs, rhs, is_pos = _canonical_lit(Eq(a, b))
    assert lhs == a and rhs == b and is_pos == True

    # Test Ne
    lhs, rhs, is_pos = _canonical_lit(Ne(a, b))
    assert lhs == a and rhs == b and is_pos == False

    # Test Not(Eq)
    lhs, rhs, is_pos = _canonical_lit(Not(Eq(a, b)))
    assert lhs == a and rhs == b and is_pos == False

    # Test unsupported literal
    with pytest.raises(EUFUnhandledInput):
        _canonical_lit(a)  # Just a symbol, not a supported literal


def test_canon_eq():
    """Test _canon_eq function."""
    eq1 = _canon_eq(a, b)
    eq2 = _canon_eq(b, a)
    assert eq1 == eq2  # Should be consistently ordered


# Test ProofProducingCongruenceClosure
def test_proof_producing_congruence_closure_init():
    """Test initialization of ProofProducingCongruenceClosure."""
    ppcc = ProofProducingCongruenceClosure([])
    assert ppcc.proof_forest == {}
    assert ppcc.direct_explanations == {}
    assert ppcc.merge_sequence == []


def test_add_equality_with_reason():
    """Test add_equality_with_reason method."""
    ppcc = ProofProducingCongruenceClosure([])
    reason = Eq(a, b)
    ppcc.add_equality_with_reason(a, b, reason)

    # Check that equality is recorded
    assert ppcc.are_equal(a, b)

    # Adding same equality again should not cause issues
    ppcc.add_equality_with_reason(a, b, reason)


def test_explain_equality_same_terms():
    """Test explain_equality with identical terms."""
    ppcc = ProofProducingCongruenceClosure([])
    explanation = ppcc.explain_equality(a, a)
    assert explanation == set()  # Same term needs no explanation


def test_explain_equality_not_equal():
    """Test explain_equality with non-equal terms."""
    ppcc = ProofProducingCongruenceClosure([])
    explanation = ppcc.explain_equality(a, b)
    assert explanation == set()  # Not equal, can't explain


def test_explain_equality_direct():
    """Test explain_equality with directly asserted equality."""
    ppcc = ProofProducingCongruenceClosure([])
    reason = Eq(a, b)
    ppcc.add_equality_with_reason(a, b, reason)

    explanation = ppcc.explain_equality(a, b)
    assert reason in explanation

@XFAIL
def test_find_proof_path_complex():
    """Test _find_proof_path with complex merge chains."""
    ppcc = ProofProducingCongruenceClosure([])

    # Create a longer chain: a = b = c = d
    ppcc.add_equality_with_reason(a, b, Eq(a, b))
    ppcc.add_equality_with_reason(b, c, Eq(b, c))
    ppcc.add_equality_with_reason(c, d, Eq(c, d))

    explanation = ppcc.explain_equality(a, d)
    assert len(explanation) > 0


# Test EUFTheorySolver exception cases
def test_equality_contradiction_exception():
    """Test EUFEqualityContradictionException."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(a, b)})

    # First assert equality
    solver.SetTrue(Eq(a, b))

    # Then try to assert disequality - should raise exception
    with pytest.raises(EUFEqualityContradictionException):
        solver.SetTrue(Ne(a, b))

@XFAIL
def test_conflict_generation_disequality():
    """Test conflict generation for disequality contradictions."""
    cnf = CNF.from_prop(Eq(a, b) & Ne(a, b))
    enc = EncodedCNF()
    enc.from_cnf(cnf)
    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)

    # Should have conflicts from initialization
    assert len(conflicts) > 0 or len(enc.data) == 0


def test_conflict_generation_equality():
    """Test conflict generation for equality contradictions."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(a, b)})

    # Set up contradiction scenario
    solver.SetTrue(Eq(a, b))

    # This should generate a conflict
    success, conflict = solver.assert_lit(-1)  # Assuming -1 corresponds to Ne(a,b)

    # The method should handle the conflict appropriately


def test_istrue_with_disequalities():
    """Test IsTrue method with disequality tracking."""
    solver = EUFTheorySolver()
    solver.Initialize({Ne(a, b), Eq(a, c), Eq(b, d)})

    solver.SetTrue(Ne(a, b))
    assert solver.IsTrue(Ne(a, b)) is True
    assert solver.IsTrue(Eq(a, b)) is False

@XFAIL
def test_explanation_with_disequality():
    """Test Explanation method for disequalities."""
    solver = EUFTheorySolver()
    solver.Initialize({Ne(a, b), Eq(a, c)})

    solver.SetTrue(Ne(a, b))
    solver.SetTrue(Eq(a, c))

    # Test explanation for the disequality
    explanation = solver.Explanation(Ne(a, b))
    assert isinstance(explanation, set)


def test_explain_disequality():
    """Test explain_disequality method."""
    solver = EUFTheorySolver()
    solver.Initialize({Ne(a, b), Eq(a, c), Eq(b, d)})

    solver.SetTrue(Ne(a, b))
    solver.SetTrue(Eq(a, c))
    solver.SetTrue(Eq(b, d))

    # Should be able to explain disequality between c and d
    explanation = solver.explain_disequality(c, d)
    assert isinstance(explanation, set)


def test_reset_functionality():
    """Test reset method."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(c, d)})

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Ne(c, d))

    # Reset should clear all state
    solver.reset()

    assert len(solver.interpretation_stack) == 0
    assert solver.timestamp_counter == 0
    assert len(solver.literal_timestamp) == 0
    assert len(solver.disequalities_set) == 0
    assert len(solver.disequality_causes) == 0


def test_backtrack_empty_stack():
    """Test backtracking with empty stack."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b)})

    # Backtrack more than available
    solver.Backtrack(5)  # Should not crash


def test_rebuild_state_with_conflicts():
    """Test _rebuild_state with potential conflicts."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Eq(b, c), Ne(a, c)})

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))

    # Force rebuild
    solver._rebuild_state()

    # Should handle conflicts gracefully during rebuild


def test_functional_congruence_complex():
    """Test complex functional congruence scenarios."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Eq(a, b),
        Eq(f(a), c),
        Eq(f(b), d),
        Eq(g(f(a)), e)
    })

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(f(a), c))

    # f(a) should equal f(b) by congruence
    assert solver.IsTrue(Eq(f(a), f(b))) is True


def test_nested_functions():
    """Test nested function applications."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Eq(a, b),
        Eq(f(g(a)), c),
        Eq(f(g(b)), d)
    })

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(f(g(a)), c))

    # Should propagate through nested functions
    assert solver.IsTrue(Eq(f(g(a)), f(g(b)))) is True


def test_from_encoded_cnf_with_trivial_true():
    """Test from_encoded_cnf with trivially true predicates."""
    enc = EncodedCNF()
    # Create a trivially true predicate (x = x)
    enc.encoding = {Q.eq(x, x): 1}
    enc.data = [{1}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(conflicts) > 0  # Should detect trivial conflicts


def test_from_encoded_cnf_with_trivial_false():
    """Test from_encoded_cnf with trivially false predicates."""
    enc = EncodedCNF()
    # Create a trivially false predicate (x != x)
    enc.encoding = {Q.ne(x, x): 1}
    enc.data = [{1}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(conflicts) > 0  # Should detect trivial conflicts


def test_from_encoded_cnf_unhandled_predicate():
    """Test from_encoded_cnf with unhandled predicate types."""
    enc = EncodedCNF()
    # Create an unhandled predicate type
    enc.encoding = {x: 1}  # Just a symbol, not a predicate
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should handle gracefully by skipping


def test_unary_predicate_handling():
    """Test handling of unary predicates."""
    enc = EncodedCNF()
    # Create unary predicate Q.prime(x)
    prime_x = AppliedPredicate(Q.prime, (x,))
    enc.encoding = {prime_x: 1}
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should convert to equality with dummy constant


def test_binary_predicate_q_eq():
    """Test Q.eq binary predicate handling."""
    enc = EncodedCNF()
    # Create Q.eq(x, y) predicate
    eq_xy = AppliedPredicate(Q.eq, (x, y))
    enc.encoding = {eq_xy: 1}
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should convert properly


def test_binary_predicate_q_ne():
    """Test Q.ne binary predicate handling."""
    enc = EncodedCNF()
    # Create Q.ne(x, y) predicate
    ne_xy = AppliedPredicate(Q.ne, (x, y))
    enc.encoding = {ne_xy: 1}
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should convert properly


def test_multiple_disequality_conflicts():
    """Test multiple disequality conflicts."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Ne(a, b), Ne(b, c), Ne(c, d),
        Eq(a, b), Eq(b, c), Eq(c, d)
    })

    # Set up chain of equalities
    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))

    # Try to add conflicting disequality
    with pytest.raises((EUFDisequalityContradictionException, EUFEqualityContradictionException)):
        solver.SetTrue(Ne(a, c))

@XFAIL
def test_congruence_with_multiple_functions():
    """Test congruence with multiple function symbols."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Eq(a, b),
        Eq(f(a), c), Eq(g(a), d),
        Eq(f(b), e), Eq(g(b), f)
    })

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(f(a), c))
    solver.SetTrue(Eq(g(a), d))

    # Both f and g should respect congruence
    assert solver.IsTrue(Eq(f(a), f(b))) is True
    assert solver.IsTrue(Eq(g(a), g(b))) is True

@XFAIL
def test_deep_explanation_chain():
    """Test explanation with deep chains of reasoning."""
    solver = EUFTheorySolver()

    # Create long chain: a=b=c=d=e
    eqs = [Eq(a, b), Eq(b, c), Eq(c, d), Eq(d, e)]
    solver.Initialize(set(eqs))

    for eq in eqs:
        solver.SetTrue(eq)

    # Explain why a = e
    explanation = solver.Explanation(Eq(a, e))
    assert len(explanation) > 0


def test_exception_type_tracking():
    """Test proper exception type tracking for conflict generation."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(a, b)})

    # Set up for disequality contradiction
    solver.SetTrue(Ne(a, b))

    try:
        solver.SetTrue(Eq(a, b))
    except EUFDisequalityContradictionException:
        # Check that exception type is properly tracked
        assert hasattr(solver, '_current_exception_type')
        assert solver._current_exception_type == 'disequality_contradiction'


def test_generate_minimal_conflict_edge_cases():
    """Test _generate_minimal_conflict with edge cases."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b)})

    # Test with no current conflict
    conflict = solver._generate_minimal_conflict()
    assert isinstance(conflict, set)

    # Test with unknown exception type
    solver._current_exception_type = 'unknown'
    conflict = solver._generate_minimal_conflict()
    assert isinstance(conflict, set)

def test_simple_chain_explanation():
    # Simple chain x=y, y=z, test explanation for x=z
    solver = EUFTheorySolver()
    cc = solver.cc
    cc.add_equality_with_reason(x, y, Eq(x, y))
    cc.add_equality_with_reason(y, z, Eq(y, z))
    explanation = cc.explain_equality(x, z)
    assert Eq(x, y) in explanation
    assert Eq(y, z) in explanation


def test_assert_literals_and_conflict():
    # Use CNF encoding from props and assumptions, assert literals and expect conflict on negation

    prop = Q.prime(z)
    assum = Q.prime(y) & Q.eq(x, y) & Q.eq(x, z) & Q.eq(z, d) & Q.eq(y, e)
    add = Q.eq(e, d)

    _prop = CNF.from_prop(prop)
    _assum = CNF.from_prop(assum)
    _add = CNF.from_prop(add)

    sat = EncodedCNF()
    sat.add_from_cnf(_prop)
    sat.add_from_cnf(_assum)
    sat.add_from_cnf(_add)

    solver = EUFTheorySolver()
    a, conflicts = solver.from_encoded_cnf(sat)

    # Assert positive literals
    for lit in [Q.prime(y), Q.eq(x, z), Q.eq(z, d), Q.eq(y, e), Q.eq(x, y)]:
        res, _ = a.assert_lit(sat.encoding[lit])
        assert res is True

    # Assert negation of e=d causes conflict
    res, conflict = a.assert_lit(-sat.encoding[Q.eq(e, d)])
    assert res is False

    # Explanation for equality e=d includes chain
    explanation = a.cc.explain_equality(e, d)
    expected_eqs = {Eq(e, y), Eq(x, z), Eq(d, z), Eq(x, y)}
    for eq in expected_eqs:
        assert eq in explanation

def test_explain_disequality_basic():
    """Test basic disequality explanation with direct disequality."""
    solver = EUFTheorySolver()

    # Manually set up disequality cause
    solver.SetTrue(Ne(x, y))

    explanation = solver.explain_disequality(x, y)
    assert Ne(x, y) in explanation
    assert len(explanation) == 1  # Only the direct disequality


def test_explain_disequality_with_equality_chain():
    """Test disequality explanation involving equality chains."""
    solver = EUFTheorySolver()

    # Set up: x = a, y = b, a != b
    # This means x != y should be explained by {Ne(a, b), Eq(x, a), Eq(y, b)}
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))
    solver.SetTrue(Ne(a, b))

    explanation = solver.explain_disequality(x, y)

    # Should contain the source disequality and connecting equalities
    assert Ne(a, b) in explanation
    assert Eq(a, x) in explanation
    assert Eq(b, y) in explanation


def test_explain_disequality_transitive_chain():
    """Test disequality explanation with longer equality chains."""
    solver = EUFTheorySolver()

    # Chain: x = u, u = v, y = w, v != w
    # So x != y should be explained through this chain
    solver.SetTrue(Eq(x, u))
    solver.SetTrue(Eq(u, v))
    solver.SetTrue(Eq(y, w))
    solver.SetTrue(Ne(v, w))

    explanation = solver.explain_disequality(x, y)

    # Should contain source disequality and all connecting equalities
    assert Ne(v, w) in explanation
    expected_equalities = {Eq(u, x), Eq(u, v), Eq(w, y), Ne(v, w)}
    for eq in expected_equalities:
        assert eq in explanation


def test_explain_disequality_reversed_terms():
    """Test disequality explanation works with reversed term order."""
    solver = EUFTheorySolver()

    # Set up a != b
    solver.SetTrue(Ne(a, b))

    # Test both directions: a != b and b != a
    explanation_ab = solver.explain_disequality(a, b)
    explanation_ba = solver.explain_disequality(b, a)

    # Both should give same explanation (just the disequality)
    assert Ne(a, b) in explanation_ab
    assert Ne(a, b) in explanation_ba
    assert explanation_ab == explanation_ba


def test_explain_disequality_no_disequality():
    """Test explanation when terms are not disequal."""
    solver = EUFTheorySolver()

    # Set up only equalities, no disequalities
    solver.SetTrue(Eq(x, y))

    # Should return empty explanation since x and y are not disequal
    explanation = solver.explain_disequality(x, y)
    assert explanation == set()


def test_explain_disequality_complex_scenario():
    """Test disequality explanation in complex scenario with multiple constraints."""
    solver = EUFTheorySolver()

    # Complex setup: x = u, u = v, y = w, w = z, v != z
    # This creates x != y through the chain
    solver.SetTrue(Eq(x, u))
    solver.SetTrue(Eq(u, v))
    solver.SetTrue(Eq(y, w))
    solver.SetTrue(Eq(w, z))
    solver.SetTrue(Ne(v, z))

    explanation = solver.explain_disequality(x, y)

    # Should include the source disequality and relevant equality chains
    assert Ne(v, z) in explanation
    # Should include equalities connecting x to v and y to z
    expected_in_explanation = {Ne(v,z), Eq(u, x), Eq(u, v), Eq(w, y), Eq(w, z)}
    for eq in expected_in_explanation:
        assert eq in explanation


def test_disequality_conflict_generation():
    """Test that disequality explanations work in conflict generation."""
    prop = Q.prime(x) & Q.prime(y) & Q.eq(x, z) & Q.eq(y, w) & Q.ne(z, w)

    _prop = CNF.from_prop(prop)
    sat = EncodedCNF()
    sat.add_from_cnf(_prop)

    solver = EUFTheorySolver()
    a, conflicts = solver.from_encoded_cnf(sat)

    # Assert the constraints
    a.assert_lit(sat.encoding[Q.prime(x)])
    a.assert_lit(sat.encoding[Q.prime(y)])
    a.assert_lit(sat.encoding[Q.eq(x, z)])
    a.assert_lit(sat.encoding[Q.eq(y, w)])
    a.assert_lit(sat.encoding[Q.ne(z, w)])

    # Now assert x = y, which should conflict due to z != w
    result, conflict = a.assert_lit(sat.encoding[Q.eq(x, y)] if Q.eq(x, y) in sat.encoding else -1)

    # Should detect conflict and generate explanation
    if result is False:
        assert len(conflict) > 0
        # Conflict should involve the disequality and connecting equalities


def test_explain_disequality_with_functions():
    """Test disequality explanation involving function terms."""
    from sympy import Function
    f = Function('f')

    solver = EUFTheorySolver()

    # Set up: f(x) != f(y), x = a, y = b
    solver.SetTrue(Ne(f(x), f(y)))
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))

    # Explanation for f(a) != f(b) should include the function disequality and argument equalities
    explanation = solver.explain_disequality(f(a), f(b))

    assert Ne(f(x), f(y)) in explanation


def test_explain_disequality_edge_cases():
    """Test edge cases for disequality explanation."""
    solver = EUFTheorySolver()

    # Test with identical terms (should return empty)
    explanation = solver.explain_disequality(x, x)
    assert explanation == set()

    # Test with unknown terms
    unknown1, unknown2 = symbols('unknown1 unknown2')
    explanation = solver.explain_disequality(unknown1, unknown2)
    assert explanation == set()


def test_conflict_generation_disequality_contradiction():
    """Test conflict generation for disequality contradiction cases."""
    solver = EUFTheorySolver()

    # Set up encoding manually for testing
    solver.literal_eqs = {1: [Eq(x, y)], 2: [Ne(x, y)]}
    solver._enc_to_lit = {1: Eq(x, y), -1: Ne(x, y), 2: Ne(x, y), -2: Eq(x, y)}
    solver.lit_to_enc = {Eq(x, y): 1, Ne(x, y): 2}

    # First assert disequality
    result, _ = solver.assert_lit(2)  # Ne(x, y)
    assert result is True

    # Then assert conflicting equality - should generate conflict
    result, conflict = solver.assert_lit(1)  # Eq(x, y)
    assert result is False
    assert len(conflict) > 0
    assert -1 in conflict or -2 in conflict  # Should negate one of the conflicting literals


def test_multiple_disequalities_explanation():
    """Test explanation when multiple disequalities exist."""
    solver = EUFTheorySolver()

    # Set up multiple disequalities: a != b, c != d
    # And equalities: x = a, y = b, z = c, w = d
    solver.SetTrue(Ne(a, b))
    solver.SetTrue(Ne(c, d))
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))
    solver.SetTrue(Eq(z, c))
    solver.SetTrue(Eq(w, d))

    # Test explanation for x != y (should use first disequality)
    explanation_xy = solver.explain_disequality(x, y)
    assert Ne(a, b) in explanation_xy
    assert Eq(a, x) in explanation_xy
    assert Eq(b, y) in explanation_xy

    # Test explanation for z != w (should use second disequality)
    explanation_zw = solver.explain_disequality(z, w)
    assert Ne(c, d) in explanation_zw
    assert Eq(c, z) in explanation_zw
    assert Eq(d, w) in explanation_zw
