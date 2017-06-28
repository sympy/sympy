from __future__ import print_function, division

from sympy import S
from sympy.core.basic import Basic
from sympy.core.function import Lambda
from sympy.core.logic import fuzzy_bool
from sympy.logic.boolalg import And
from sympy.sets.sets import (Set, Interval, Intersection, EmptySet, Union,
                             FiniteSet)
from sympy.utilities.iterables import sift


class ConditionSet(Set):
    """
    Set of elements which satisfies a given condition.

    {x | condition(x) is True for x in S}

    Examples
    ========

    >>> from sympy import Symbol, S, ConditionSet, Lambda, pi, Eq, sin, Interval
    >>> x = Symbol('x')
    >>> sin_sols = ConditionSet(x, Eq(sin(x), 0), Interval(0, 2*pi))
    >>> 2*pi in sin_sols
    True
    >>> pi/2 in sin_sols
    False
    >>> 3*pi in sin_sols
    False
    >>> 5 in ConditionSet(x, x**2 > 4, S.Reals)
    True
    """
    def __new__(cls, sym, condition, base_set):
        if condition == S.false:
            return S.EmptySet
        if condition == S.true:
            return base_set
        if isinstance(base_set, EmptySet):
            return base_set
        if isinstance(base_set, FiniteSet):
            sifted = sift(base_set, lambda _: fuzzy_bool(condition.subs(sym, _)))
            if sifted[None]:
                return Union(FiniteSet(*sifted[True]),
                             Basic.__new__(cls, sym, condition, FiniteSet(*sifted[None])))
            else:
                return FiniteSet(*sifted[True])
        return Basic.__new__(cls, sym, condition, base_set)

    sym = property(lambda self: self.args[0])
    condition = property(lambda self: self.args[1])
    base_set = property(lambda self: self.args[2])

    def _intersect(self, other):
        if not isinstance(other, ConditionSet):
            return ConditionSet(self.sym, self.condition,
                                Intersection(self.base_set, other))

    def contains(self, other):
        return And(Lambda(self.sym, self.condition)(other), self.base_set.contains(other))
