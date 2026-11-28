"""
rename this to test_assumptions.py when the old assumptions system is deleted
"""
from __future__ import annotations
from sympy.abc import x, y
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.ask import Q, _ask_single_fact
from sympy.printing import pretty
from sympy.assumptions.cnf import CNF


def test_equal():
    """Test for equality"""
    assert Q.positive(x) == Q.positive(x)
    assert Q.positive(x) != ~Q.positive(x)
    assert ~Q.positive(x) == ~Q.positive(x)


def test_pretty():
    assert pretty(Q.positive(x)) == "Q.positive(x)"
    assert pretty(
        {Q.positive, Q.integer}) == "{Q.integer, Q.positive}"


def test_global():
    """Test for global assumptions"""
    global_assumptions.add(x > 0)
    assert (x > 0) in global_assumptions
    global_assumptions.remove(x > 0)
    assert not (x > 0) in global_assumptions
    # same with multiple of assumptions
    global_assumptions.add(x > 0, y > 0)
    assert (x > 0) in global_assumptions
    assert (y > 0) in global_assumptions
    global_assumptions.clear()
    assert not (x > 0) in global_assumptions
    assert not (y > 0) in global_assumptions


def test_ask_single_fact():
    assert _ask_single_fact(Q.real, CNF()) is None
    assert _ask_single_fact(Q.even, CNF.from_prop(Q.zero)) is True
    assert _ask_single_fact(Q.even, CNF.from_prop(Q.odd)) is False
    assert _ask_single_fact(Q.even, CNF.from_prop(Q.real)) is None
    assert _ask_single_fact(Q.integer, CNF.from_prop(Q.even | Q.odd)) is None
    assert _ask_single_fact(Q.integer, CNF.from_prop(Q.prime)) is True
    assert _ask_single_fact(Q.prime,   CNF.from_prop(Q.composite)) is False
    assert _ask_single_fact(Q.zero, CNF.from_prop(~Q.even)) is False

    assert _ask_single_fact(Q.zero, CNF.from_prop(~Q.even & Q.real)) is False

def test_encoded_cnf_representation_and_add_clause():
    from sympy.assumptions.cnf import EncodedCNF
    from sympy.core.singleton import S
    from sympy.testing.pytest import raises

    cnf = EncodedCNF()
    assert cnf.data == []

    cnf.add_prop(S.true)
    assert cnf.data == []

    # Tautologies and duplicate literals are normalized/dropped
    cnf.add_clause([1, -1])
    assert cnf.data == []
    cnf.add_clause([1, 1, -2])
    cnf.add_clause([-2, 1])
    assert cnf.data == [{1, -2}]

    # Invalid literals raise ValueError
    raises(ValueError, lambda: cnf.add_clause([0]))
    raises(ValueError, lambda: cnf.add_clause([True]))
    raises(ValueError, lambda: cnf.add_clause(["x"]))

    # Empty clause collapses CNF to [set()] and blocks further additions
    cnf.add_clause([])
    assert cnf.data == [set()]
    cnf.add_clause([2])
    assert cnf.data == [set()]

    # S.false encodes to [set()], never [{0}]
    false_cnf = EncodedCNF()
    false_cnf.add_prop(S.false)
    assert false_cnf.data == [set()]
    assert {0} not in false_cnf.data
