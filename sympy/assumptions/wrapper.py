"""
Functions and wrapper object to call assumption property and predicate
query with same syntax.

In SymPy, there are two assumption systems. Old assumption system is
defined in sympy/core/assumptions, and it can be accessed by attribute
such as ``x.is_even``. New assumption system is definded in
sympy/assumptions, and it can be accessed by predicates such as
``Q.even(x)``.

Old assumption is fast, while new assumptions can freely take local facts.
In general, old assumption is used in evaluation method and new assumption
is used in refinement method.

In most cases, both evaluation and refinement follow the same process, and
the only difference is which assumption system is used. This module provides
``is_[...]()`` functions and ``AssumptionsWrapper()`` class which allows
using two systems with same syntax so that parallel code implementation can be
avoided.

Examples
========

For multiple use, use ``AssumptionsWrapper()``.

>>> from sympy import Q, Symbol
>>> from sympy.assumptions.wrapper import AssumptionsWrapper
>>> x = Symbol('x')
>>> _x = AssumptionsWrapper(x, Q.even(x))
>>> _x.is_integer
True
>>> _x.is_odd
False

For single use, use ``is_[...]()`` functions.

>>> from sympy.assumptions.wrapper import is_infinite
>>> a = Symbol('a')
>>> print(is_infinite(a))
None
>>> is_infinite(a, Q.finite(a))
False

"""

from sympy.assumptions import ask, Q
from sympy.core.assumptions import (_assume_defined, as_property,
    ManagedProperties)
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify

class AssumptionsWrapperMeta(ManagedProperties):
    """
    Metaclass to give _eval_is_[...] attributes to AssumptionsWrapper
    """
    def __init__(cls, *args, **kws):
        for fact in _assume_defined:
            pname = "_eval_%s" % as_property(fact)
            setattr(cls, pname, make_eval_method(fact))
        super().__init__(cls, *args, **kws)


def make_eval_method(fact):
    def getit(self):
        try:
            pred = getattr(Q, fact)
            ret = ask(pred(self.expr), self.assumptions)
            return ret
        except AttributeError:
            return None
    return getit


# we subclass Basic to use the fact deduction and caching
class AssumptionsWrapper(Basic, metaclass=AssumptionsWrapperMeta):
    """
    Wrapper over ``Basic`` instances to call predicate query by
    ``.is_[...]`` property

    Parameters
    ==========

    expr : Basic

    assumptions : Boolean, optional

    Examples
    ========

    >>> from sympy import Q, Symbol
    >>> from sympy.assumptions.wrapper import AssumptionsWrapper
    >>> x = Symbol('x', even=True)
    >>> AssumptionsWrapper(x).is_integer
    True
    >>> y = Symbol('y')
    >>> AssumptionsWrapper(y, Q.even(y)).is_integer
    True

    With ``AssumptionsWrapper``, both evaluation and refinement can be supported
    by single implementation.

    >>> from sympy import Function
    >>> class MyAbs(Function):
    ...     @classmethod
    ...     def eval(cls, x, assumptions=True):
    ...         _x = AssumptionsWrapper(x, assumptions)
    ...         if _x.is_nonnegative:
    ...             return x
    ...         if _x.is_negative:
    ...             return -x
    ...     def _eval_refine(self, assumptions):
    ...         return MyAbs.eval(self.args[0], assumptions)
    >>> MyAbs(x)
    MyAbs(x)
    >>> MyAbs(x).refine(Q.positive(x))
    x
    >>> MyAbs(Symbol('y', negative=True))
    -y

    """
    def __new__(cls, expr, assumptions=None):
        if assumptions is None:
            return expr
        obj = super().__new__(cls, expr, _sympify(assumptions))
        obj.expr = expr
        obj.assumptions = assumptions
        return obj


# one shot functions which are faster than AssumptionsWrapper

def is_infinite(obj, assumptions=None):
    if assumptions is None:
        return obj.is_infinite
    return ask(Q.infinite(obj), assumptions)


def is_extended_real(obj, assumptions=None):
    if assumptions is None:
        return obj.is_extended_real
    return ask(Q.extended_real(obj), assumptions)


def is_extended_nonnegative(obj, assumptions=None):
    if assumptions is None:
        return obj.is_extended_nonnegative
    return ask(Q.extended_nonnegative(obj), assumptions)
