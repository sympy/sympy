from sympy.sets import ConditionSet
from sympy import (Symbol, Eq, S, sin, pi, Lambda)

x = Symbol('x')


def test_CondSet():
    sin_sols = ConditionSet(Lambda(x, Eq(sin(x), 0)), S.Reals)
    assert 2*pi in sin_sols
