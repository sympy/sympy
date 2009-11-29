"""rename this to test_assumptions.py when the old assumptions system is deleted"""
from sympy.core import symbols
from sympy.assumptions import Assume, global_assumptions
from sympy.assumptions.assume import eliminate_assume
from sympy.printing import pretty

def test_assume():
    x = symbols('x')
    assump = Assume(x, 'integer')
    assert assump.expr == x
    assert assump.key == 'integer'
    assert assump.value == True

def test_False():
    """Test Assume object with False keys"""
    x = symbols('x')
    assump = Assume(x, 'integer', False)
    assert assump.expr == x
    assert assump.key == 'integer'
    assert assump.value == False

def test_equal():
    """Test for equality"""
    x = symbols('x')
    assert Assume(x, 'positive', True)  == Assume(x, 'positive', True)
    assert Assume(x, 'positive', True)  != Assume(x, 'positive', False)
    assert Assume(x, 'positive', False) == Assume(x, 'positive', False)

def test_pretty():
    x = symbols('x')
    assert pretty(Assume(x, 'positive', True)) == "Assume(x, 'positive', True)"

def test_eliminate_assumptions():
    a, b, x, y = symbols('abxy')
    assert eliminate_assume(Assume(x, 'a', True))  == a
    assert eliminate_assume(Assume(x, 'a', True), symbol=x)  == a
    assert eliminate_assume(Assume(x, 'a', True), symbol=y)  == None
    assert eliminate_assume(Assume(x, 'a', False)) == ~a
    assert eliminate_assume(Assume(x, 'a', False), symbol=y) == None
    assert eliminate_assume(Assume(x, 'a', True) | Assume(x, 'b')) == a | b
    assert eliminate_assume(Assume(x, 'a', True) | Assume(x, 'b', False)) == a | ~b

def test_global():
    """Test for global assumptions"""
    x, y = symbols('x y')
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
