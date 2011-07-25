"""rename this to test_assumptions.py when the old assumptions system is deleted"""
from sympy.core import symbols
from sympy.assumptions import AppliedPredicate, global_assumptions, Predicate
from sympy.assumptions.ask import _extract_facts
from sympy.printing import pretty
from sympy.assumptions.ask import Q
from sympy.logic.boolalg import Or
from sympy.utilities.pytest import XFAIL

def test_equal():
    """Test for equality"""
    x = symbols('x')
    assert Q.positive(x)  == Q.positive(x)
    assert Q.positive(x)  != ~Q.positive(x)
    assert ~Q.positive(x) == ~Q.positive(x)

def test_pretty():
    x = symbols('x')
    assert pretty(Q.positive(x)) == "Q.positive(x)"

def test_extract_facts():
    a, b = symbols('a b', cls=Predicate)
    x, y = symbols('x y')
    assert _extract_facts(a(x), x)  == a
    assert _extract_facts(a(x), y)  == None
    assert _extract_facts(~a(x), x) == ~a
    assert _extract_facts(~a(x), y) == None
    assert _extract_facts(a(x) | b(x), x) == a | b
    assert _extract_facts(a(x) | ~b(x), x) == a | ~b

def test_global():
    """Test for global assumptions"""
    x, y = symbols('x,y')
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

def test_composite_predicates():
    x = symbols('x')
    pred = Q.integer | ~Q.positive
    assert type(pred(x)) is Or
    assert pred(x) == Q.integer(x) | ~Q.positive(x)
