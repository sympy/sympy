from sympy.sets import (ConditionSet, Intersection)
from sympy import (Symbol, Eq, S, Abs, sin, pi, Lambda, Interval)

x = Symbol('x')


def test_CondSet():
    sin_sols_principal = ConditionSet(Lambda(x, Eq(sin(x), 0)),
                                      Interval(0, 2*pi, False, True))
    assert pi in sin_sols_principal
    assert pi/2 not in sin_sols_principal
    assert 3*pi not in sin_sols_principal
    assert 5 in ConditionSet(Lambda(x, x**2 > 4), S.Reals)
    assert 1 not in ConditionSet(Lambda(x, x**2 > 4), S.Reals)


def test_CondSet_intersect():
    input_conditionset = ConditionSet(Lambda(x, x**2 > 4), Interval(1, 4, False, False))
    other_domain = Interval(0, 3, False, False)
    output_conditionset = ConditionSet(Lambda(x, x**2 > 4), Interval(1, 3, False, False))
    assert Intersection(input_conditionset, other_domain) == output_conditionset


def test_issue_9849():
    assert ConditionSet(Lambda(x, Eq(x, x)), S.Naturals) == S.Naturals
    assert ConditionSet(Lambda(x, Eq(Abs(sin(x)), -1)), S.Naturals) == S.EmptySet
