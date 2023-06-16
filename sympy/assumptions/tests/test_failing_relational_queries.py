from sympy.testing.pytest import XFAIL
from sympy.assumptions.ask import ask, Q
from sympy.logic.boolalg import Equivalent, Implies, And, Or
from sympy.abc import x, y, z



@XFAIL
def test_inequalities():
    # test basic defintions
    assert ask(Equivalent(Q.gt(x, 0), Q.positive(x))) is True
    assert ask(Equivalent(Q.eq(x, 0), Q.zero(x))) is True
    assert ask(Equivalent(Q.lt(x, 0), Q.negative(x))) is True
    assert ask(Equivalent(Q.ge(x, 0), Or(Q.positive(x), Q.zero(x)))) is True
    assert ask(Equivalent(Q.le(x, 0), Or(Q.negative(x), Q.zero(x)))) is True
    assert ask(Implies(And(Q.gt(x, y), Q.positive(y)), Q.positive(x))) is True
    assert ask(Implies(And(Q.lt(x, y), Q.negative(y)), Q.negative(x))) is True

    # test more complex problems
    assert ask(x > z, x > y & y > z)


@XFAIL
def test_equality():
    # test substitution property of equality
    assert ask(Q.prime(x), Q.eq(x, y) & Q.prime(y)) is True
    assert ask(Q.real(x), Q.eq(x, y) & Q.real(y)) is True
    assert ask(Q.imaginary(x), Q.eq(x, y) & Q.imaginary(y)) is True

@XFAIL
def test_mixed():
    # test assumptions that mix inequalities and assumptions that the sat solver handles
    assert ask(x>0, Q.positive(x) & Q.prime(y))




