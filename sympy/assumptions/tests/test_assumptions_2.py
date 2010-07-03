"""rename this to test_assumptions.py when the old assumptions system is deleted"""
from sympy.core import symbols
from sympy.assumptions import Assume, global_assumptions, Predicate
from sympy.assumptions.assume import eliminate_assume
from sympy.printing import pretty
from sympy.assumptions.ask import Q
from sympy.utilities.pytest import XFAIL

def test_assume():
    x = symbols('x')
    assump = Assume(x, 'integer')
    assert assump.expr == x
    assert assump.key == Q.integer

def test_Predicate_wraps_Assume():
    x = symbols('x')
    integer = Predicate('integer')
    assump = integer(x)
    assert (assump.expr, assump.key) == (x, integer)
    assump = Assume(x, integer)
    assert (assump.expr, assump.key) == (x, integer)

def test_False():
    """Test Assume object with False keys"""
    x = symbols('x')
    assump = Assume(x, 'integer', False)
    assert assump == ~Assume(x, 'integer')


def test_equal():
    """Test for equality"""
    x = symbols('x')
    assert Assume(x, 'positive', True)  == Assume(x, 'positive', True)
    assert Assume(x, 'positive', True)  != Assume(x, 'positive', False)
    assert Assume(x, 'positive', False) == Assume(x, 'positive', False)

@XFAIL #TODO: handle printing
def test_pretty():
    x = symbols('x')
    assert pretty(Assume(x, 'positive')) == "Assume(x, 'positive')"

def test_eliminate_assumptions():
    a, b = map(Predicate, symbols('a,b'))
    x, y = symbols('x,y')
    assert eliminate_assume(Assume(x, a))  == a
    assert eliminate_assume(Assume(x, a), symbol=x)  == a
    assert eliminate_assume(Assume(x, a), symbol=y)  == None
    assert eliminate_assume(Assume(x, a, False)) == ~a
    assert eliminate_assume(Assume(x, a), symbol=y) == None
    assert eliminate_assume(Assume(x, a) | Assume(x, b)) == a | b
    assert eliminate_assume(Assume(x, a) | Assume(x, b, False)) == a | ~b

def test_global():
    """Test for global assumptions"""
    x, y = symbols('x,y')
    global_assumptions.add(Assume(x>0))
    assert Assume(x>0) in global_assumptions
    global_assumptions.remove(Assume(x>0))
    assert not Assume(x>0) in global_assumptions
    # same with multiple of assumptions
    global_assumptions.add(Assume(x>0), Assume(y>0))
    assert Assume(x>0) in global_assumptions
    assert Assume(y>0) in global_assumptions
    global_assumptions.clear()
    assert not Assume(x>0) in global_assumptions
    assert not Assume(y>0) in global_assumptions
