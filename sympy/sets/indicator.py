from sympy.core.symbol import Symbol
from sympy.core.function import Function
from sympy.ntheory.residue_ntheory import mobius

from .sets import Set


class indicator(Function):
    r"""
    Indicator function (also called the Characteristic function) maps an element, x, to the range {0, 1} for set A

    It is defined as follows:
        1) `1` if x exists in A
        2) `0` if x does not exist in A

    Parameters
    ==========

    x : the element
    a : the set

    Examples
    ========

    >>> from sympy.sets import Union, Interval, FiniteSet, indicator
    >>> set = Union(Interval(0, 1), Interval(2, 3))
    >>> indicator(0.5, set)
    1
    >>> indicator(1.5, set)
    0
    >>> indicator(2.5, set)
    1
    >>> set = FiniteSet(1, 2, 3)
    >>> indicator(0, set)
    0
    >>> indicator(1, set)
    1


    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Indicator_function
    .. [2] https://mathworld.wolfram.com/CharacteristicFunction.html

    """

    @classmethod
    def eval(cls, x, a: Set):
        if not isinstance(x, Symbol):
            if type(a) != Set:
                return 1 if a.contains(x) else 0

    def inverse(self):
        """
        Returns the Mobius function, the inverse of the indicator function. [1]
        """
        return mobius
