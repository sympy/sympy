from sympy.sets import (ConditionSet, Intersection, FiniteSet, EmptySet, Union)
from sympy import (Symbol, Eq, S, Abs, sin, pi, Lambda, Interval, And, Mod)

x = Symbol('x')


def test_CondSet():
    sin_sols_principal = ConditionSet(x, Eq(sin(x), 0),
                                      Interval(0, 2*pi, False, True))
    assert pi in sin_sols_principal
    assert pi/2 not in sin_sols_principal
    assert 3*pi not in sin_sols_principal
    assert 5 in ConditionSet(x, x**2 > 4, S.Reals)
    assert 1 not in ConditionSet(x, x**2 > 4, S.Reals)


def test_CondSet_intersect():
    input_conditionset = ConditionSet(x, x**2 > 4, Interval(1, 4, False, False))
    other_domain = Interval(0, 3, False, False)
    output_conditionset = ConditionSet(x, x**2 > 4, Interval(1, 3, False, False))
    assert Intersection(input_conditionset, other_domain) == output_conditionset


def test_issue_9849():
    assert ConditionSet(x, Eq(x, x), S.Naturals) == S.Naturals
    assert ConditionSet(x, Eq(Abs(sin(x)), -1), S.Naturals) == S.EmptySet

def test_simplified_FiniteSet_in_CondSet():
    assert ConditionSet(x, And(x < 1, x > -3), FiniteSet(0, 1, 2)) == FiniteSet(0)
    assert ConditionSet(x, x < 0, FiniteSet(0, 1, 2)) == EmptySet()
    assert ConditionSet(x, And(x < -3), EmptySet()) == EmptySet()
    y = Symbol('y')
    assert (ConditionSet(x, And(x > 0), FiniteSet(-1, 0, 1, y)) ==
        Union(FiniteSet(1), ConditionSet(x, And(x > 0), FiniteSet(y))))
    assert (ConditionSet(x, Eq(Mod(x, 3), 1), FiniteSet(1, 4, 2, y)) ==
        Union(FiniteSet(1, 4), ConditionSet(x, Eq(Mod(x, 3), 1), FiniteSet(y))))
