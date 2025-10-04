"""
Functions and wrapper object to call assumption property and predicate
query with same syntax.

In SymPy, there are two assumption systems. Old assumption system is
defined in sympy/core/assumptions, and it can be accessed by attribute
such as ``x.is_even``. New assumption system is defined in
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
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.singleton import S

import itertools


def make_eval_method(fact):
    def getit(self):
        pred = getattr(Q, fact)
        ret = ask(pred(self.expr), self.assumptions)
        return ret
    return getit


# we subclass Basic to use the fact deduction and caching
class AssumptionsWrapper(Basic):
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

    _eval_is_algebraic = make_eval_method("algebraic")
    _eval_is_antihermitian = make_eval_method("antihermitian")
    _eval_is_commutative = make_eval_method("commutative")
    _eval_is_complex = make_eval_method("complex")
    _eval_is_composite = make_eval_method("composite")
    _eval_is_even = make_eval_method("even")
    _eval_is_extended_negative = make_eval_method("extended_negative")
    _eval_is_extended_nonnegative = make_eval_method("extended_nonnegative")
    _eval_is_extended_nonpositive = make_eval_method("extended_nonpositive")
    _eval_is_extended_nonzero = make_eval_method("extended_nonzero")
    _eval_is_extended_positive = make_eval_method("extended_positive")
    _eval_is_extended_real = make_eval_method("extended_real")
    _eval_is_finite = make_eval_method("finite")
    _eval_is_hermitian = make_eval_method("hermitian")
    _eval_is_imaginary = make_eval_method("imaginary")
    _eval_is_infinite = make_eval_method("infinite")
    _eval_is_integer = make_eval_method("integer")
    _eval_is_irrational = make_eval_method("irrational")
    _eval_is_negative = make_eval_method("negative")
    _eval_is_noninteger = make_eval_method("noninteger")
    _eval_is_nonnegative = make_eval_method("nonnegative")
    _eval_is_nonpositive = make_eval_method("nonpositive")
    _eval_is_nonzero = make_eval_method("nonzero")
    _eval_is_odd = make_eval_method("odd")
    _eval_is_polar = make_eval_method("polar")
    _eval_is_positive = make_eval_method("positive")
    _eval_is_prime = make_eval_method("prime")
    _eval_is_rational = make_eval_method("rational")
    _eval_is_real = make_eval_method("real")
    _eval_is_transcendental = make_eval_method("transcendental")
    _eval_is_zero = make_eval_method("zero")


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


def assumption_domain(symbol, assumptions=set(), _exclude=frozenset()):
    from .assume import global_assumptions

    from sympy.core.exprtools import _monotonic_sign, _eps
    from sympy.core.symbol import Symbol, Dummy
    from sympy.sets import FiniteSet, Interval, imageset

    if not isinstance(symbol, Symbol):
        raise ValueError("Input must be a symbol, use calculus.util.function_domain() for expressions")

    if symbol.is_infinite and symbol.is_extended_positive:
        return FiniteSet(S.Infinity)
    elif symbol.is_infinite and symbol.is_extended_negative:
        return FiniteSet(S.NegativeInfinity)
    elif symbol.is_zero:
        return FiniteSet(S.Zero)

    domain = S.UniversalSet

    for assumption in itertools.chain(global_assumptions, assumptions):
        if domain == S.EmptySet:
            break

        if assumption in _exclude or not assumption.has(symbol):
            continue

        domain = domain.intersect(assumption.as_set(symbol, _exclude=_exclude | {symbol}))

    if symbol.is_negative is False:
        result = _monotonic_sign(symbol)
        if result is not None:
            if result in (_eps, -_eps):
                domain = domain.intersect(Interval(S.Zero, S.Infinity, S.true, S.false))
            else:
                domain = domain.intersect(Interval(result, S.Infinity, S.false, S.false))
    elif symbol.is_positive is False:
        result = _monotonic_sign(symbol)
        if result is not None:
            if result in (_eps, -_eps):
                domain = domain.intersect(Interval(S.NegativeInfinity, S.Zero, S.false, S.true))
            else:
                domain = domain.intersect(Interval(S.NegativeInfinity, result, S.false, S.false))

    if symbol.is_even:
        n = Dummy('n')
        domain = domain.intersect(imageset(n, n * 2, S.Integers))
    elif symbol.is_odd:
        n = Dummy('n')
        domain = domain.intersect(imageset(n, n * 2 + 1, S.Integers))

    # TODO: add case for is_prime

    return domain
