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

# For now, we join all assumptions together, so it makes the most sense for
# things assumed in old and new to be conjuncted together. This can lead to
# situations where ask gives False simply because the assumptions themselves
# are inconsistent (which is actually better than the current system, see
# issue 3917).

def test_assuming_old_and_new_asking_new1():
    x = Symbol('x', positive=True)
    with assuming(Q.negative(x)):
        assert ask(Q.negative(x)) is False

    assert ask(Q.negative(x)) is False
    assert ask(Q.negative(x), Q.negative(x)) is False

def test_assuming_old_and_new_asking_old1():
    x = Symbol('x', positive=True)
    with assuming(Q.negative(x)):
        assert x.is_negative is False

# It is different if we ask on something that was assumed in the old system or
# if we ask on something that was assumed in the new one.
def test_assuming_old_and_new_asking_new2():
    x = Symbol('x', positive=True)
    with assuming(Q.negative(x)):
        assert ask(Q.positive(x)) is False

    assert ask(Q.positive(x), Q.negative(x)) is False

def test_assuming_old_and_new_asking_old2():
    x = Symbol('x', positive=True)
    with assuming(Q.negative(x)):
        assert x.is_positive is False
