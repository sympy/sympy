"""This tests sympy/core/basic.py with (ideally) no reference to subclasses
of Basic or Atom."""

from sympy.core.basic import Basic, Atom, S

from sympy.utilities.pytest import raises


b1 = Basic(); b2 = Basic(b1); b3 = Basic(b2)
b21 = Basic(b2, b1)

def test_structure():
    assert b21.args == (b2, b1)
    assert tuple(b21.iter_basic_args()) == b21.args
    assert b21.func(*b21.args) == b21
    assert b21.new(*b21.args) == b21
    assert bool(b1)

def test_equality():
    instances = [b1, b2, b3, b21, Basic(b1,b1,b1), Basic]
    for i, b_i in enumerate(instances):
        for j, b_j in enumerate(instances):
            assert (b_i == b_j) == (i == j)
            assert (b_i != b_j) == (i != j)

    assert Basic() != []
    assert not(Basic() == [])
    assert Basic() != 0
    assert not(Basic() == 0)

def test_matches_basic():
    instances = [Basic(b1,b1,b2), Basic(b1,b2,b1), Basic(b2, b1, b1),
                    Basic(b1, b2), Basic(b2, b1), b2, b1]
    for i, b_i in enumerate(instances):
        for j, b_j in enumerate(instances):
            if i ==j:
                assert b_i.matches(b_j) == {}
            else:
                assert b_i.matches(b_j) is None
    assert b1.match(b1) == {}

def test_has():
    assert b21.has(b1)
    assert b21.has(b3, b1)
    assert b21.has(Basic)
    assert not b1.has(b21, b3)
    raises(TypeError, "b21.has()")

def test_subs():
    assert b21.subs(b2, b1) == Basic(b1, b1)
    assert b21.subs(b2, b21) == Basic(b21, b1)
    assert b3.subs(b2, b1) == b2

    assert b21.subs([(b2, b1), (b1, b2)]) == Basic(b2, b2)

    assert b21.subs({b1: b2, b2: b1}) == Basic(b2, b2)

    raises(TypeError, "b21.subs('bad arg')")
    raises(TypeError, "b21.subs(b1, b2, b3)")

def test_atoms():
    assert b21.atoms() == set()

def test_doit():
    assert b21.doit() == b21
    assert b21.doit(deep=False) == b21

def test_S():
    assert repr(S) == 'S'
