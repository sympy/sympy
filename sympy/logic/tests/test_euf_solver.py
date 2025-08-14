import pytest
from sympy import symbols, Function, Eq, Not, Unequality
from sympy.logic.algorithms.euf_solver import EUFTheorySolver
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic import boolalg

f, g = symbols('f g', cls=Function)
a, b, c, d, e, x, y, z, u = symbols('a b c d e x y z u')


def test_initialize_and_istrue():
    solver = EUFTheorySolver()
    eqs = {Eq(a, b), Eq(b, c), Eq(f(a), f(b))}
    solver.Initialize(eqs)

    # Initially false
    for lit in eqs:
        assert solver.IsTrue(lit) is False

    # Assert a=b
    solver.SetTrue(Eq(a, b))
    assert solver.IsTrue(Eq(a, b)) is True
    assert solver.IsTrue(Eq(b, a)) is True
    # others still undecided
    assert solver.IsTrue(Eq(b, c)) is False


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
    with pytest.raises(ValueError):
        solver.SetTrue(Eq(a, b))


def test_backtrack_removes_last():
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Eq(b, c)})
    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))
    assert solver.IsTrue(Eq(a, c))
    solver.Backtrack(1)
    assert solver.IsTrue(Eq(a, c)) is False
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


def test_initialize_type_error():
    solver = EUFTheorySolver()
    with pytest.raises(TypeError):
        solver.Initialize({a})  # not Eq or Not(Eq) or Unequality


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

def test_basic_initialization_and_encoding():
    solver, enc, conflicts = make_solver_from_props(Eq(x, y))
    assert conflicts == []
    # Mapping enc_id -> literal should be non-empty
    assert len(solver._enc_to_lit) > 0
    # Check that x != x is trivially false conflict
    cnf = CNF.from_prop(Eq(x, x))
    enc2 = EncodedCNF(); enc2.from_cnf(cnf)
    solver2, conflicts2 = EUFTheorySolver.from_encoded_cnf(enc2)
    assert conflicts2  # should detect trivial equality

def test_assert_lit_positive_and_check_sat():
    solver, enc, conflicts = make_solver_from_props(Eq(x, y))
    assert not conflicts
    # Get the encoded literal ID (positive int)
    lit_id = next(iter(solver._enc_to_lit))
    # Nothing asserted yet → IsTrue() is False
    assert solver.IsTrue(solver._enc_to_lit[lit_id]) is False
    # Assert as positive literal
    assert solver.assert_lit(lit_id) is None
    # Should now be seen as True
    assert solver.IsTrue(solver._enc_to_lit[lit_id]) is True
    # check() should still be SAT
    sat, info = solver.check()
    assert sat
    assert isinstance(info, dict)

def test_order_independence_of_assertions():
    # x=y, y=z — test that order doesn't matter
    solver, enc, conflicts = make_solver_from_props(Eq(x, y), Eq(y, z))
    ids = list(solver._enc_to_lit)
    solver.assert_lit(ids[1])
    solver.assert_lit(ids[0])
    assert solver.IsTrue(Eq(x, z)) is True
    # Reset and reverse
    solver.Backtrack(len(solver.interpretation_stack))
    solver.assert_lit(ids[0])
    solver.assert_lit(ids[1])
    assert solver.IsTrue(Eq(x, z)) is True

def test_simple_equality_chain():
    """
    Test EUF: x = y, y = z  => SAT, and x = z must hold.
    """
    cnf = CNF().from_prop(Eq(x, y) & Eq(y, z))
    enc = EncodedCNF(); enc.from_cnf(cnf)
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)

    # Assert all clauses/literals
    for lit_id in enc.encoding.values():
        assert euf.assert_lit(lit_id) is None

    # Check satisfiability
    is_sat, _ = euf.check()
    assert is_sat is True
    # Derived fact
    assert euf.IsTrue(Eq(x, z)) is True

def test_backtrack_recovery():
    """
    EUF: assert eq, then backtrack and verify it’s gone.
    """
    cnf = CNF().from_prop(Eq(x, y))
    enc = EncodedCNF(); enc.from_cnf(cnf)
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    lit_id = next(iter(enc.encoding.values()))

    euf.assert_lit(lit_id)
    assert euf.IsTrue(Eq(x, y)) is True

    euf.Backtrack(1)  # undo the only decision
    assert euf.IsTrue(Eq(x, y)) is False
    sat, _ = euf.check()
    assert sat is True
