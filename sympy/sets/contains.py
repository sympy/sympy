from __future__ import print_function, division

from sympy.core import Basic
from sympy.core.relational import Eq, Ne
from sympy.logic.boolalg import BooleanFunction


class Contains(BooleanFunction):
    """
    Asserts that x is an element of the set S

    Examples
    ========

    >>> from sympy import Symbol, Integer, S
    >>> from sympy.sets.contains import Contains
    >>> Contains(Integer(2), S.Integers)
    True
    >>> Contains(Integer(-2), S.Naturals)
    False
    >>> i = Symbol('i', integer=True)
    >>> Contains(i, S.Naturals)
    Contains(i, S.Naturals)

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Element_%28mathematics%29
    """
    @classmethod
    def eval(cls, x, S):
        from sympy.sets.sets import Set

        if not isinstance(x, Basic):
            raise TypeError
        if not isinstance(S, Set):
            raise TypeError

        ret = S.contains(x)
        if not isinstance(ret, Contains):
            return ret

    @property
    def binary_symbols(self):
        return set().union(*[i.binary_symbols
            for i in self.args[1].args
            if i.is_Boolean or i.is_Symbol or
            isinstance(i, (Eq, Ne))])

    def as_set(self):
        return self
