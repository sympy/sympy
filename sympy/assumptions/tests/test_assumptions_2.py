"""
rename this to test_assumptions.py when the old assumptions system is deleted
"""
from sympy.abc import x, y
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.ask import Q
from sympy.printing import pretty


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
    global_assumptions.add(Q.is_true(x > 0))
    assert Q.is_true(x > 0) in global_assumptions
    global_assumptions.remove(Q.is_true(x > 0))
    assert not Q.is_true(x > 0) in global_assumptions
    # same with multiple of assumptions
    global_assumptions.add(Q.is_true(x > 0), Q.is_true(y > 0))
    assert Q.is_true(x > 0) in global_assumptions
    assert Q.is_true(y > 0) in global_assumptions
    global_assumptions.clear()
    assert not Q.is_true(x > 0) in global_assumptions
    assert not Q.is_true(y > 0) in global_assumptions
