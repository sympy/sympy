from sympy.sets import ConditionSet
from sympy import (Symbol, Eq, S, sin, pi, Lambda, Interval)

x = Symbol('x')


def test_CondSet():
    sin_sols_principal = ConditionSet(Lambda(x, Eq(sin(x), 0)),
                                      Interval(0, 2*pi, False, True))
    assert pi in sin_sols_principal
    assert pi/2 not in sin_sols_principal
    assert 3*pi not in sin_sols_principal
    assert 5 in ConditionSet(Lambda(x, x**2 > 4), S.Reals)
    assert 1 not in ConditionSet(Lambda(x, x**2 > 4), S.Reals)
