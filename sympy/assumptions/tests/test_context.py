from __future__ import annotations
from sympy.assumptions import ask, Q
from sympy.assumptions.assume import assuming, global_assumptions
from sympy.abc import x, y

def test_assuming():
    with assuming(Q.integer(x)):
        assert ask(Q.integer(x))
    assert not ask(Q.integer(x))

def test_assuming_nested():
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))
    with assuming(Q.integer(x)):
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
        with assuming(Q.integer(y)):
            assert ask(Q.integer(x))
            assert ask(Q.integer(y))
        assert ask(Q.integer(x))
        assert not ask(Q.integer(y))
    assert not ask(Q.integer(x))
    assert not ask(Q.integer(y))

def test_finally():
    try:
        with assuming(Q.integer(x)):
            1/0
    except ZeroDivisionError:
        pass
    assert not ask(Q.integer(x))

def test_remove_safe():
    global_assumptions.add(Q.integer(x))
    with assuming():
        assert ask(Q.integer(x))
        global_assumptions.remove(Q.integer(x))
        assert not ask(Q.integer(x))
    assert ask(Q.integer(x))
    global_assumptions.clear() # for the benefit of other tests

def test_encoded_cnf_invariants():
    from sympy.assumptions.cnf import EncodedCNF
    from sympy.testing.pytest import raises

    enc = EncodedCNF()
    assert enc.data == []

    raises(ValueError, lambda: enc.add_clause({0}))
    raises(ValueError, lambda: enc.add_clause({True}))
    raises(ValueError, lambda: enc.add_clause({"x"}))

    # Tautologies are omitted
    enc.add_clause({1, -1})
    assert enc.data == []

    # Valid clauses are deduplicated
    enc.add_clause({1, 2})
    enc.add_clause({1, 2})
    assert enc.data == [{1, 2}]

    # Empty clause collapses CNF to [set()] (False)
    enc.add_clause(set())
    assert enc.data == [set()]

    # Further clauses after contradiction do not alter [set()]
    enc.add_clause({1})
    assert enc.data == [set()]


def test_encoded_cnf_boolean_constants():
    from sympy.assumptions.cnf import CNF, EncodedCNF
    from sympy.core.singleton import S

    enc = EncodedCNF()
    enc.add_from_cnf(CNF.from_prop(S.true))
    assert enc.data == []

    enc.add_from_cnf(CNF.from_prop(S.false))
    assert enc.data == [set()]
