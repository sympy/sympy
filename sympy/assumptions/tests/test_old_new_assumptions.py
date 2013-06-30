from __future__ import with_statement
from sympy import Symbol, assuming, Q, ask, S

def test_assuming_new_asking_old():
    x = Symbol('x')
    with assuming(Q.positive(x)):
        assert x.is_positive

def test_assuming_old_asking_new():
    x = Symbol('x', positive=True)
    assert ask(Q.positive(x))

def test_assuming_old_asking_old():
    x = Symbol('x', positive=True)
    assert x.is_positive

def test_assuming_new_asking_new():
    x = Symbol('x')
    with assuming(Q.positive(x)):
        assert ask(Q.positive(x))

def test_defined_asking_new():
    assert ask(Q.positive(2))

def test_defined_asking_old():
    assert S(2).is_positive
