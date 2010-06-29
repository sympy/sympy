from sympy import symbols, Integer
from sympy.tensor import Indexed, Idx
from sympy.utilities.pytest import raises
from sympy import oo

def test_Idx_construction():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i) != Idx(i, 1)
    assert Idx(i, a) == Idx(i, (0, a - 1))
    assert Idx(i, oo) == Idx(i, (0, oo))

    x = symbols('x')
    raises(TypeError, "Idx(x)")
    raises(TypeError, "Idx(0.5)")
    raises(TypeError, "Idx(i, x)")
    raises(TypeError, "Idx(i, 0.5)")
    raises(TypeError, "Idx(i, (x, 5))")
    raises(TypeError, "Idx(i, (2, x))")
    raises(TypeError, "Idx(i, (2, 3.5))")

def test_Idx_properties():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i).is_integer

def test_Idx_bounds():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i).lower == None
    assert Idx(i).upper == None
    assert Idx(i, a).lower == 0
    assert Idx(i, a).upper == a - 1
    assert Idx(i, 5).lower == 0
    assert Idx(i, 5).upper == 4
    assert Idx(i, oo).lower == 0
    assert Idx(i, oo).upper == oo
    assert Idx(i, (a, b)).lower == a
    assert Idx(i, (a, b)).upper == b
    assert Idx(i, (1, 5)).lower == 1
    assert Idx(i, (1, 5)).upper == 5
    assert Idx(i, (-oo, oo)).lower == -oo
    assert Idx(i, (-oo, oo)).upper == oo

def test_Idx_fixed_bounds():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(2).lower == None
    assert Idx(2).upper == None
    assert Idx(2, a).lower == 0
    assert Idx(2, a).upper == a - 1
    assert Idx(2, 5).lower == 0
    assert Idx(2, 5).upper == 4
    assert Idx(2, oo).lower == 0
    assert Idx(2, oo).upper == oo
    assert Idx(2, (a, b)).lower == a
    assert Idx(2, (a, b)).upper == b
    assert Idx(2, (1, 5)).lower == 1
    assert Idx(2, (1, 5)).upper == 5
    assert Idx(2, (-oo, oo)).lower == -oo
    assert Idx(2, (-oo, oo)).upper == oo

def test_Idx_func_args():
    i, a, b = symbols('i a b', integer=True)
    ii = Idx(i)
    assert ii.func(*ii.args) == ii
    ii = Idx(i, a)
    assert ii.func(*ii.args) == ii
    ii = Idx(i, (a, b))
    assert ii.func(*ii.args) == ii

def test_Idx_subs():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i, a).subs(a, b) == Idx(i, b)
    assert Idx(i, a).subs(i, b) == Idx(b, a)

    assert Idx(i).subs(i,2) == Idx(2)
    assert Idx(i, a).subs(a, 2) == Idx(i, 2)
    assert Idx(i, (a, b)).subs(i, 2) == Idx(2, (a, b))
