from sympy.sets.condset import CondSet
from sympy import (Symbol, Eq, S, sin, pi, Lambda)

x = Symbol('x')


def test_CondSet():
    sin_sols = CondSet(Lambda(x, Eq(sin(x), 0)), S.Reals)
    assert 2*pi in sin_sols
