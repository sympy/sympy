"""
Singularities
=============

This module implements algorithms for finding singularities for a function
and identifying types of functions.

The differential calculus methods in this module include methods to identify
the following function types in the given ``Interval``:
- Increasing
- Strictly Increasing
- Decreasing
- Strictly Decreasing
- Monotonic

"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify
from sympy.functions.elementary.exponential import log
from sympy.functions.elementary.trigonometric import sec, csc, cot, tan, cos
from sympy.functions.elementary.hyperbolic import (
    sech, csch, coth, tanh, cosh, asech, acsch, atanh, acoth)
from sympy.utilities.misc import filldedent


if TYPE_CHECKING:
    from typing import Callable
    from sympy.core.expr import Expr
    from sympy.sets.sets import Set
    from sympy.core.basic import Basic
    from sympy.logic.boolalg import Boolean


def singularities(
    expression: Expr | complex,
    symbol: Symbol | Basic,
    domain: Set | None = None,
) -> Set:
    """
    Find singularities of a given function.
    """
    from sympy.solvers.solveset import solveset
    from sympy.sets.sets import Interval

    expression = sympify(expression)

    if domain is None:
        domain = S.Reals if symbol.is_real else S.Complexes
    try:
        sings = S.EmptySet
        e = expression.rewrite([sec, csc, cot, tan], cos)
        e = e.rewrite([sech, csch, coth, tanh], cosh)
        for i in e.atoms(Pow):
            if i.base == S.Zero:
                sing_interval = solveset(i.exp <= 0, symbol, domain)
                if isinstance(sing_interval, Interval):
                    sing_interval = Interval.open(sing_interval.inf, sing_interval.sup)
                sings += sing_interval
            if i.exp.is_infinite:
                raise NotImplementedError
            if i.exp.is_negative:
                sings += solveset(i.base, symbol, domain)
        for i in expression.atoms(log, asech, acsch):
            sings += solveset(i.args[0], symbol, domain)
        for i in expression.atoms(atanh, acoth):
            sings += solveset(i.args[0] - 1, symbol, domain)
            sings += solveset(i.args[0] + 1, symbol, domain)
        return sings
    except NotImplementedError:
        raise NotImplementedError(filldedent('''
            Methods for determining the singularities
            of this function have not been developed.'''))


###########################################################################
#                     DIFFERENTIAL CALCULUS METHODS                       #
###########################################################################


def monotonicity_helper(
    expression: Expr | complex,
    predicate: Callable[[Expr], Boolean],
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Helper function for checking function monotonicity.
    """
    from sympy.solvers.solveset import solveset

    expression = sympify(expression)
    free = expression.free_symbols

    if symbol is None:
        if len(free) > 1:
            raise NotImplementedError(
                'The function has not yet been implemented'
                ' for all multivariate expressions.'
            )

    variable = symbol or (free.pop() if free else Symbol('x'))

    try:
        sings = singularities(expression, variable, interval)
    except NotImplementedError:
        pass
    else:
        if interval.is_subset(S.Reals):
            interior_sings = interval.interior.intersection(sings)
            if interior_sings != S.EmptySet:
                return False

    derivative = expression.diff(variable)
    predicate_interval = solveset(predicate(derivative), variable, S.Reals)
    return interval.is_subset(predicate_interval)


def is_increasing(
    expression: Expr | complex,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """Return whether the function is increasing in the given interval."""
    return monotonicity_helper(expression, lambda x: x >= 0, interval, symbol)


def is_strictly_increasing(
    expression: Expr | complex,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """Return whether the function is strictly increasing in the given interval."""
    return monotonicity_helper(expression, lambda x: x > 0, interval, symbol)


def is_decreasing(
    expression: Expr | complex,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """Return whether the function is decreasing in the given interval."""
    return monotonicity_helper(expression, lambda x: x <= 0, interval, symbol)


def is_strictly_decreasing(
    expression: Expr | complex,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """Return whether the function is strictly decreasing in the given interval."""
    return monotonicity_helper(expression, lambda x: x < 0, interval, symbol)


def is_monotonic(
    expression: Expr | complex,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Return whether the function is monotonic in the given interval.
    A function is monotonic if it is either non-decreasing or non-increasing
    across the entire specified interval.
    """
    expression = sympify(expression)
    free = expression.free_symbols

    if symbol is None and len(free) > 1:
        raise NotImplementedError(
            'is_monotonic has not yet been implemented'
            ' for all multivariate expressions.'
        )

    return (
        is_increasing(expression, interval, symbol) or
        is_decreasing(expression, interval, symbol)
    )
