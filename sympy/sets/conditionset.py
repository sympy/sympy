from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.logic.boolalg import And, Or, Not, true, false
from sympy.sets.sets import (Set, Interval, Intersection, EmptySet, Union,
                             FiniteSet)
from sympy.core.singleton import Singleton, S
from sympy.core.sympify import _sympify
from sympy.core.decorators import deprecated
from sympy.core.function import Lambda


class ConditionSet(Set):
    """
    Set of elements which satisfies a given condition.

    {x | condition(x) is True for x in S}

    Examples
    ========

    >>> from sympy import Symbol, S, CondSet, Lambda, pi, Eq, sin

    >>> x = Symbol('x')
    >>> sin_sols = CondSet(Lambda(x, Eq(sin(x), 0)), S.Reals)
    >>> 2*pi in sin_sols
    True
    """
    def __new__(cls, lamda, base_set):
        return Basic.__new__(cls, lamda, base_set)

    condition = property(lambda self: self.args[0])
    base_set = property(lambda self: self.args[1])

    def contains(self, other):
        return And(self.condition(other), self.base_set.contains(other))
