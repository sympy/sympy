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


def test_issue_27834():
    """Regression test for sympy/sympy#27834.

    `ask(Q.zero(expr), Q.eq(var, value))` should give the same result whether
    the variable is declared as ``real=True`` at construction or whether the
    realness is added later via ``Q.real(var)`` (or omitted entirely, since
    Q.eq(v, 1) implies v == 1 which is real).
    """
    from sympy import ask, symbols

    # Case 1: real=True at construction -- this worked before the fix.
    x = symbols('x', real=True)
    assert ask(Q.zero(x - 1), Q.eq(x, 1)) is True

    # Case 2: Q.real added in the assumption clause -- this is what the
    # original bug report complained about.
    y = symbols('y')
    assert ask(Q.zero(y - 1), Q.eq(y, 1) & Q.real(y)) is True

    # Case 3: Q.eq on its own (realness is implied by the equality).
    z = symbols('z')
    assert ask(Q.zero(z - 1), Q.eq(z, 1)) is True

    # Case 4: order of Q.real and Q.eq should not matter.
    w = symbols('w')
    assert ask(Q.zero(w - 1), Q.real(w) & Q.eq(w, 1)) is True

    # Sanity: ask() should still return False for expressions that are NOT
    # forced to zero by the assumption.
    assert ask(Q.zero(y - 2), Q.eq(y, 1)) is False
    assert ask(Q.zero(y - 2), Q.eq(y, 1) & Q.real(y)) is False
