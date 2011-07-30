"""This tests sympy/core/basic.py with (ideally) no reference to subclasses
of Basic or Atom."""

from sympy.core.basic import Basic, Atom
from sympy.core.singleton import S, Singleton

from sympy.utilities.pytest import raises


b1 = Basic(); b2 = Basic(b1); b3 = Basic(b2)
b21 = Basic(b2, b1)

def test_structure():
    assert b21.args == (b2, b1)
    assert tuple(b21.iter_basic_args()) == b21.args
    assert b21.func(*b21.args) == b21
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
    assert not b21.has()

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

def test_free_symbols_empty():
    assert b21.free_symbols == set()

def test_doit():
    assert b21.doit() == b21
    assert b21.doit(deep=False) == b21

def test_S():
    assert repr(S) == 'S'

def test_Singleton():
    global instanciated
    instanciated = 0
    class MySingleton(Basic):
        __metaclass__ = Singleton

        def __new__(cls):
            global instanciated
            instanciated += 1
            return Basic.__new__(cls)

    assert instanciated == 1
    assert MySingleton() is not Basic()
    assert MySingleton() is MySingleton()
    assert S.MySingleton is MySingleton()
    assert instanciated == 1

    class MySingleton_sub(MySingleton):
        pass
    assert instanciated == 2
    assert MySingleton_sub() is not MySingleton()
    assert MySingleton_sub() is MySingleton_sub()
