"""
sympy/logic/tests/test_euf_theory.py

Tests for the EUF (Equality with Uninterpreted Functions) theory solver.
"""
from sympy.assumptions.ask import Q
from sympy.logic.algorithms.euf_theory import EUFSolver, EUFUnhandledInput
from sympy import symbols, Symbol
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask_generated import get_known_facts_dict
from sympy.testing.pytest import raises

x, y, z = symbols("x y z")

def test_euf_basic_predicates():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.even, True)
    assert solver.query(Q.even, x) is True
    assert solver.query(Q.odd, x) is False
    assert solver.query(Q.integer, x) is True
    solver.assign_fact(solver.get_var(y), Q.odd, True)
    assert solver.query(Q.odd, y) is True
    assert solver.query(Q.even, y) is False


def test_euf_composite_and_prime():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.composite, True)
    assert solver.query(Q.prime, x) is False
    assert solver.query(Q.integer, x) is True
    assert solver.query(Q.positive, x) is True
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.prime, True)
    assert solver.query(Q.composite, x) is False
    assert solver.query(Q.integer, x) is True
    assert solver.query(Q.positive, x) is True


def test_euf_rational_irrational_real():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.rational, True)
    assert solver.query(Q.real, x) is True
    assert solver.query(Q.irrational, x) is False
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.irrational, True)
    assert solver.query(Q.real, x) is True
    assert solver.query(Q.rational, x) is False


def test_euf_matrix_properties():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.orthogonal, True)
    assert solver.query(Q.unitary, x) is True
    assert solver.query(Q.positive_definite, x) is True
    assert solver.query(Q.invertible, x) is True
    assert solver.query(Q.fullrank, x) is True
    assert solver.query(Q.square, x) is True
    assert solver.query(Q.singular, x) is False


def test_euf_equality_and_disequality():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.even, True)
    solver.add_equality(x, y)
    assert solver.query(Q.even, y) is True
    solver.assign_fact(solver.get_var(z), Q.odd, True)
    solver.add_disequality(x, z)
    assert solver.check_consistency() is True
    raises(ValueError, lambda: solver.add_equality(y, z))


def test_euf_reflexivity_and_transitivity():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.even, True)
    solver.add_equality(x, y)
    solver.add_equality(y, z)
    assert solver.query(Q.even, z) is True


def test_euf_matrix_exclusivity():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.invertible, True)
    assert solver.query(Q.singular, x) is False
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.singular, True)
    assert solver.query(Q.invertible, x) is False


def test_euf_negative_propagation():
    # Negative propagation: assigning False propagates negations
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.composite, False)
    # Q.composite False should propagate Q.prime True (by exclusivity)
    assert solver.query(Q.prime, x) is True
    # Q.composite False should propagate Q.integer False, Q.positive False, etc.
    assert solver.query(Q.integer, x) is False
    assert solver.query(Q.positive, x) is False


def test_euf_all_fact_propagation():
    facts = get_known_facts_dict()
    for pred, (implied, rejected) in facts.items():
        s = EUFSolver()
        s.assign_fact(s.get_var(x), pred, True)
        for p in implied:
            assert s.query(p, x) is True, f"{pred} should imply {p}"
        for p in rejected:
            assert s.query(p, x) is False, f"{pred} should reject {p}"


def test_euf_all_fact_negative_propagation():
    facts = get_known_facts_dict()
    for pred, (implied, rejected) in facts.items():
        s = EUFSolver()
        s.assign_fact(s.get_var(x), pred, False)
        for p in implied:
            assert s.query(p, x) is False, f"~{pred} should reject {p}"
        for p in rejected:
            assert s.query(p, x) is True, f"~{pred} should imply {p}"


def test_euf_large_equivalence_class():
    # Test propagation through a chain of equalities
    vars = symbols('a b c d e f')
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(vars[0]), Q.even, True)
    for i in range(len(vars) - 1):
        solver.add_equality(vars[i], vars[i+1])
    for v in vars:
        assert solver.query(Q.even, v) is True
        assert solver.query(Q.integer, v) is True


def test_euf_disequality_conflict():
    solver = EUFSolver()
    solver.add_equality(x, y)
    solver.add_disequality(x, y)
    assert solver.check_consistency() is False


def test_euf_unknown_predicate():
    solver = EUFSolver()
    assert solver.query(Q.transcendental, x) is None


def test_euf_multiple_predicates():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.rational, True)
    solver.assign_fact(solver.get_var(x), Q.positive, True)
    assert solver.query(Q.real, x) is True
    assert solver.query(Q.integer, x) is None


def test_euf_matrix_chain():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.orthogonal, True)
    # Orthogonal implies unitary, invertible, fullrank, square, positive_definite, etc.
    assert solver.query(Q.unitary, x) is True
    assert solver.query(Q.invertible, x) is True
    assert solver.query(Q.fullrank, x) is True
    assert solver.query(Q.square, x) is True
    assert solver.query(Q.positive_definite, x) is True
    assert solver.query(Q.singular, x) is False


def test_euf_equivalence():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.irrational, True)
    assert solver.query(Q.real, x) is True
    assert solver.query(Q.rational, x) is False


def test_euf_reflexivity():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.positive, True)
    assert solver.query(Q.positive, x) is True


def test_euf_exclusivity_matrix():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.invertible, True)
    assert solver.query(Q.singular, x) is False
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.singular, True)
    assert solver.query(Q.invertible, x) is False


def test_euf_multiple_disequalities():
    a, b, c = symbols('a b c')
    solver = EUFSolver()
    solver.add_equality(a, b)
    solver.add_disequality(b, c)
    solver.assign_fact(solver.get_var(a), Q.even, True)
    solver.assign_fact(solver.get_var(c), Q.odd, True)
    assert solver.check_consistency() is True
    raises(ValueError, lambda: solver.add_equality(a, c))



def test_euf_chain_implications():
    solver = EUFSolver()
    solver.assign_fact(solver.get_var(x), Q.orthogonal, True)
    # Should propagate to all implied matrix properties
    assert solver.query(Q.unitary, x) is True
    assert solver.query(Q.invertible, x) is True
    assert solver.query(Q.fullrank, x) is True
    assert solver.query(Q.square, x) is True
    assert solver.query(Q.positive_definite, x) is True
    assert solver.query(Q.singular, x) is False


def test_euf_full_fact_table():
    # Exhaustively test all predicates in the known facts dict
    facts = get_known_facts_dict()
    for pred in facts:
        s = EUFSolver()
        s.assign_fact(s.get_var(x), pred, True)
        implied, rejected = facts[pred]
        for p in implied:
            assert s.query(p, x) is True
        for p in rejected:
            assert s.query(p, x) is False


def test_euf_negative_full_fact_table():
    # Exhaustively test all predicates in the known facts dict (negative)
    facts = get_known_facts_dict()
    x = Symbol('x')
    for pred in facts:
        s = EUFSolver()
        s.assign_fact(s.get_var(x), pred, False)
        implied, rejected = facts[pred]
        for p in implied:
            assert s.query(p, x) is False
        for p in rejected:
            assert s.query(p, x) is True


def build_solver_from_assumptions(assumptions, context=None):
    cnf = CNF.from_prop(assumptions)
    if context:
        cnf.extend(context)
    enc_cnf = EncodedCNF()
    enc_cnf.from_cnf(cnf)
    return EUFSolver.from_encoded_cnf(enc_cnf)


def test_encoded_cnf_composite_and_prime():
    solver = build_solver_from_assumptions(Q.composite(y) & Q.eq(x, y))
    assert solver.query(Q.prime, x) is False
    assert solver.query(Q.integer, x) is True
    assert solver.query(Q.positive, x) is True

    solver = build_solver_from_assumptions(Q.prime(y) & Q.eq(x, y))
    assert solver.query(Q.composite, x) is False
    assert solver.query(Q.integer, x) is True
    assert solver.query(Q.positive, x) is True


def test_encoded_cnf_rational_irrational_real():
    solver = build_solver_from_assumptions(Q.rational(y) & Q.eq(x, y))
    assert solver.query(Q.real, x) is True
    assert solver.query(Q.irrational, x) is False

    solver = build_solver_from_assumptions(Q.irrational(y) & Q.eq(x, y))
    assert solver.query(Q.real, x) is True
    assert solver.query(Q.rational, x) is False


def test_encoded_cnf_matrix_properties():
    solver = build_solver_from_assumptions(Q.orthogonal(y) & Q.eq(x, y))
    assert solver.query(Q.unitary, x) is True
    assert solver.query(Q.positive_definite, x) is True
    assert solver.query(Q.invertible, x) is True
    assert solver.query(Q.fullrank, x) is True
    assert solver.query(Q.square, x) is True
    assert solver.query(Q.singular, x) is False


def test_encoded_cnf_transitive_equality():
    solver = build_solver_from_assumptions(Q.even(x) & Q.eq(x, y) & Q.eq(y, z))
    assert solver.query(Q.even, z) is True
    assert solver.query(Q.integer, z) is True


def test_encoded_cnf_ambiguous_cases():
    solver = build_solver_from_assumptions(Q.prime(y) & Q.eq(x, y))
    # x could be 2 (even+prime) or 5 (odd+prime), so even(x) is unknown
    assert solver.query(Q.even, x) is None
    # No info about x
    solver = build_solver_from_assumptions(True)
    assert solver.query(Q.positive, x) is None

def test_encoded_cnf_large_equivalence_class():
    vars = symbols('a b c d e f')
    solver = build_solver_from_assumptions(Q.even(vars[0]) & Q.eq(vars[0], vars[1]) &
                                           Q.eq(vars[1], vars[2]) & Q.eq(vars[2], vars[3]) &
                                           Q.eq(vars[3], vars[4]) & Q.eq(vars[4], vars[5]))
    for v in vars:
        assert solver.query(Q.even, v) is True
        assert solver.query(Q.integer, v) is True
