from __future__ import print_function, division

from sympy import S
from sympy.core.basic import Basic
from sympy.logic.boolalg import And
from sympy.sets.sets import (Set, Interval, Intersection, EmptySet, Union,
                             FiniteSet)


class ConditionSet(Set):
    """
    Set of elements which satisfies a given condition.

    {x | condition(x) is True for x in S}

    Examples
    ========

    >>> from sympy import Symbol, S, ConditionSet, Lambda, pi, Eq, sin, Interval
    >>> x = Symbol('x')
    >>> sin_sols = ConditionSet(Lambda(x, Eq(sin(x), 0)), Interval(0, 2*pi))
    >>> 2*pi in sin_sols
    True
    >>> pi/2 in sin_sols
    False
    >>> 3*pi in sin_sols
    False
    >>> 5 in ConditionSet(Lambda(x, x**2 > 4), S.Reals)
    True
    """
    def __new__(cls, condition, base_set):
        if condition.args[1] is S.false:
            return S.EmptySet
        if condition.args[1] is S.true:
            return base_set
        return Basic.__new__(cls, condition, base_set)

    condition = property(lambda self: self.args[0])
    base_set = property(lambda self: self.args[1])

    def _intersect(self, other):
        if not isinstance(other, ConditionSet):
            return ConditionSet(self.condition,
                                Intersection(self.base_set, other))

    def contains(self, other):
        return And(self.condition(other), self.base_set.contains(other))
