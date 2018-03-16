from __future__ import print_function, division

from sympy import S
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.function import Lambda
from sympy.core.logic import fuzzy_bool
from sympy.core.symbol import Symbol, Dummy
from sympy.logic.boolalg import And, as_Boolean
from sympy.sets.sets import (Set, Interval, Intersection, EmptySet, Union,
                             FiniteSet)
from sympy.utilities.iterables import sift
from sympy.multipledispatch import dispatch


class ConditionSet(Set):
    """
    Set of elements which satisfies a given condition.

    {x | condition(x) is True for x in S}

    Examples
    ========

    >>> from sympy import Symbol, S, ConditionSet, pi, Eq, sin, Interval
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
        condition = as_Boolean(condition)
        if isinstance(base_set, set):
            base_set = FiniteSet(*base_set)
        elif not isinstance(base_set, Set):
            raise TypeError('expecting set for base_set')
        if condition == S.false:
            return S.EmptySet
        if condition == S.true:
            return base_set
        if isinstance(base_set, EmptySet):
            return base_set
        if isinstance(base_set, FiniteSet):
            sifted = sift(
                base_set, lambda _: fuzzy_bool(
                    condition.subs(sym, _)))
            if sifted[None]:
                return Union(FiniteSet(*sifted[True]),
                    Basic.__new__(cls, sym, condition,
                    FiniteSet(*sifted[None])))
            else:
                return FiniteSet(*sifted[True])
        if isinstance(base_set, cls):
            s, c, base_set = base_set.args
            if sym == s:
                condition = And(condition, c)
            elif sym not in c.free_symbols:
                condition = And(condition, c.xreplace({s: sym}))
            elif s not in condition.free_symbols:
                condition = And(condition.xreplace({sym: s}), c)
                sym = s
            else:
                # user will have to use cls.sym to get symbol
                dum = Symbol('lambda')
                if dum in condition.free_symbols or \
                        dum in c.free_symbols:
                    dum = Dummy(str(dum))
                condition = And(
                    condition.xreplace({sym: dum}),
                    c.xreplace({s: dum}))
                sym = dum
        if sym in base_set.free_symbols or \
                not isinstance(sym, Symbol):
            s = Symbol('lambda')
            if s in base_set.free_symbols:
                s = Dummy('lambda')
            condition = condition.xreplace({sym: s})
            sym = s
        return Basic.__new__(cls, sym, condition, base_set)

    sym = property(lambda self: self.args[0])
    condition = property(lambda self: self.args[1])
    base_set = property(lambda self: self.args[2])

    @property
    def free_symbols(self):
        s, c, b = self.args
        return (c.free_symbols - s.free_symbols) | b.free_symbols

    def contains(self, other):
        return And(Lambda(self.sym, self.condition)(
            other), self.base_set.contains(other))

    def _eval_subs(self, old, new):
        if old == self.sym:
            if new not in self.free_symbols:
                if isinstance(new, Symbol):
                    return self.func(*[i.subs(old, new) for i in self.args])
            return self.func(self.sym, self.condition, self.base_set.subs(old, new))
        return self.func(*([self.sym] + [i.subs(old, new) for i in self.args[1:]]))

    def dummy_eq(self, other, symbol=None):
        if not isinstance(other, self.func):
            return False
        if symbol:
            raise ValueError('symbol arg not supported for ConditionSet')
        o = other.func(self.sym,
            other.condition.subs(other.sym, self.sym),
            other.base_set)
        return self == o
