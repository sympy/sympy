from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.sets.sets import (Set, Interval, Intersection, EmptySet, Union,
                             FiniteSet)
from sympy.core.singleton import Singleton, S
from sympy.core.sympify import _sympify
from sympy.core.decorators import deprecated
from sympy.core.function import Lambda


class CondSet(Set):
    """
    Set of elements which satisfies a given condition.

    {x | cond(x) is True for x in S}

    Examples
    ========

    >>> from sympy import Symbol, S, CondSet, FiniteSet, Lambda, pi

    >>> x = Symbol('x')
    >>> sin_sols = CondSet(Lambda(x, Eq(sin(x), 0)), S.Reals)
    >>> 2*pi in sin_sols
    True
    """
    def __new__(cls, lamda, base_set):
        return Basic.__new__(cls, lamda, base_set)

    lamda = property(lambda self: self.args[0])
    base_set = property(lambda self: self.args[1])

    def _is_multivariate(self):
        return len(self.lamda.variables) > 1

    def _contains(self, other):
        # XXX: probably we should check if self.cond is returning only true or
        # false
        return self.cond(other)
