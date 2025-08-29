import pytest
from sympy.testing.pytest import XFAIL
from sympy import symbols, Function, Eq, Not, Unequality
from sympy.logic.algorithms.euf_theory_solver import *
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic import boolalg
from sympy.assumptions.ask import Q

f, g = symbols('f g', cls=Function)
a, b, c, d, e, x, y, z, u = symbols('a b c d e x y z u')


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
