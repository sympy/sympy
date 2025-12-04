from sympy.core import S
from sympy.core.sympify import sympify
from sympy.core.parameters import global_parameters
from sympy.logic.boolalg import Boolean
from sympy.utilities.misc import func_name
from .sets import Set


class Subset(Boolean):
    """
    Asserts that the set T is part (a subset) of the set S.

    Examples
    ========

    >>> from sympy import S, Subset, Set
    >>> Subset(S.Naturals, S.Integers)
    True
    >>> Subset(S.Reals, S.Naturals)
    False
    >>> i = Set()                           #TODO: make a SetSymbol class (like MatSymbol) and make the other set classes compatible with it.
    >>> Subset(i, S.Naturals)
    Subset(Set(), Naturals)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Subset
    """
    def __new__(cls, t, s, evaluate=None):
        t = sympify(t)
        s = sympify(s)

        if evaluate is None:
            evaluate = global_parameters.evaluate

        if not isinstance(t, Set):
            raise TypeError('expecting Set, not %s' % func_name(t))

        if not isinstance(s, Set):
            raise TypeError('expecting Set, not %s' % func_name(s))

        if evaluate:
            # is_subset can return symbolic booleans that would be returned by
            # s.is_subset(t) but here for Subset(t, s) we only evaluate to
            # true, false or return the unevaluated Contains.
            result = t.is_subset(s)
            result = sympify(result)                # apparently, the is_subset method may return python native bools...
            if isinstance(result, Boolean):
                if result in (S.true, S.false):
                    return result
            elif result is not None:
                raise TypeError("is_subset() should return Boolean or None")

        return super().__new__(cls, t, s)

    # @property
    # def binary_symbols(self):
    #     return set().union(*[i.binary_symbols
    #         for i in self.args[1].args
    #         if i.is_Boolean or i.is_Symbol or
    #         isinstance(i, (Eq, Ne))])

    # def as_set(self):
    #     return self.args[1]
