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
    expression: Expr,
    symbol: Symbol | Basic,
    domain: Set | None = None,
) -> set[Symbol]:
    """
    Find singularities of a given function.

    Parameters
    ==========

    expression : Expr
        The target function in which singularities need to be found.
    symbol : Symbol
        The symbol over the values of which the singularity in
        expression in being searched for.

    Returns
    =======

    Set
        A set of values for ``symbol`` for which ``expression`` has a
        singularity. An ``EmptySet`` is returned if ``expression`` has no
        singularities for any given value of ``Symbol``. An ``Interval`` is returned
        if ``expression`` has uncountably many singularities.

    Raises
    ======

    NotImplementedError
        Methods for determining the singularities of this function have
        not been developed.

    Notes
    =====

    This function does not find non-isolated singularities
    nor does it find branch points of the expression.

    Currently supported functions are:
        - univariate continuous (real or complex) functions

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Mathematical_singularity

    Examples
    ========

    >>> from sympy import singularities, Symbol, log
    >>> x = Symbol('x', real=True)
    >>> y = Symbol('y', real=False)
    >>> singularities(x**2 + x + 1, x)
    EmptySet
    >>> singularities(1/(x + 1), x)
    {-1}
    >>> singularities(1/(y**2 + 1), y)
    {-I, I}
    >>> singularities(1/(y**3 + 1), y)
    {-1, 1/2 - sqrt(3)*I/2, 1/2 + sqrt(3)*I/2}
    >>> singularities(log(x), x)
    {0}
    >>> singularities(0**x, x)
    Interval.open(-oo, 0)

    """
    from sympy.solvers.solveset import solveset
    from sympy.sets.sets import Interval

    if domain is None:
        domain = S.Reals if symbol.is_real else S.Complexes
    try:
        sings = S.EmptySet
        e = expression.rewrite([sec, csc, cot, tan], cos)
        e = e.rewrite([sech, csch, coth, tanh], cosh)
        for i in e.atoms(Pow):
            if i.base == S.Zero:
                sing_interval = solveset(i.exp <= 0, symbol, domain)
                # Since in Sympy we assume that 0**0 = 1,
                # we can't have a singularity at i.exp == 0 and therefore
                # we must return an open interval.
                sing_interval = Interval.open(sing_interval.inf, sing_interval.sup)
                sings += sing_interval
            if i.exp.is_infinite:
                raise NotImplementedError
            if i.exp.is_negative:
                # XXX: exponent of varying sign not handled
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
#                      DIFFERENTIAL CALCULUS METHODS                      #
###########################################################################


def monotonicity_helper(
    expression: Expr,
    predicate: Callable[[Expr], Boolean],
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Helper function for functions checking function monotonicity.

    Parameters
    ==========

    expression : Expr
        The target function which is being checked
    predicate : function
        The property being tested for. The function takes in an integer
        and returns a boolean. The integer input is the derivative and
        the boolean result should be true if the property is being held,
        and false otherwise.
    interval : Set, optional
        The range of values in which we are testing, defaults to all reals.
    symbol : Symbol, optional
        The symbol present in expression which gets varied over the given range.

    It returns a boolean indicating whether the interval in which
    the function's derivative satisfies given predicate is a superset
    of the given interval.

    Returns
    =======

    Boolean
        True if ``predicate`` is true for all the derivatives when ``symbol``
        is varied in ``range``, False otherwise.

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
        if interval.is_subset(S.Reals):
            interior_sings = interval.interior.intersection(sings)
            if interior_sings != S.EmptySet:
                return False
    except (NotImplementedError, AttributeError):
        pass

    derivative = expression.diff(variable)
    predicate_interval = solveset(predicate(derivative), variable, S.Reals)
    return interval.is_subset(predicate_interval)


def is_increasing(
    expression: Expr,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Return whether the function is increasing in the given interval.

    Parameters
    ==========

    expression : Expr
        The target function which is being checked.
    interval : Set, optional
        The range of values in which we are testing (defaults to set of
        all real numbers).
    symbol : Symbol, optional
        The symbol present in expression which gets varied over the given range.

    Returns
    =======

    Boolean
        True if ``expression`` is increasing (either strictly increasing or
        constant) in the given ``interval``, False otherwise.

    Examples
    ========

    >>> from sympy import is_increasing
    >>> from sympy.abc import x, y
    >>> from sympy import S, Interval, oo
    >>> is_increasing(x**3 - 3*x**2 + 4*x, S.Reals)
    True
    >>> is_increasing(-x**2, Interval(-oo, 0))
    True
    >>> is_increasing(-x**2, Interval(0, oo))
    False
    >>> is_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval(-2, 3))
    False
    >>> is_increasing(x**2 + y, Interval(1, 2), x)
    True

    """
    return monotonicity_helper(expression, lambda x: x >= 0, interval, symbol)


def is_strictly_increasing(
    expression: Expr,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Return whether the function is strictly increasing in the given interval.

    Parameters
    ==========

    expression : Expr
        The target function which is being checked.
    interval : Set, optional
        The range of values in which we are testing (defaults to set of
        all real numbers).
    symbol : Symbol, optional
        The symbol present in expression which gets varied over the given range.

    Returns
    =======

    Boolean
        True if ``expression`` is strictly increasing in the given ``interval``,
        False otherwise.

    Examples
    ========

    >>> from sympy import is_strictly_increasing
    >>> from sympy.abc import x, y
    >>> from sympy import Interval, oo
    >>> is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.Ropen(-oo, -2))
    True
    >>> is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.Lopen(3, oo))
    True
    >>> is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.open(-2, 3))
    False
    >>> is_strictly_increasing(-x**2, Interval(0, oo))
    False
    >>> is_strictly_increasing(-x**2 + y, Interval(-oo, 0), x)
    False

    """
    return monotonicity_helper(expression, lambda x: x > 0, interval, symbol)


def is_decreasing(
    expression: Expr,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Return whether the function is decreasing in the given interval.

    Parameters
    ==========

    expression : Expr
        The target function which is being checked.
    interval : Set, optional
        The range of values in which we are testing (defaults to set of
        all real numbers).
    symbol : Symbol, optional
        The symbol present in expression which gets varied over the given range.

    Returns
    =======

    Boolean
        True if ``expression`` is decreasing (either strictly decreasing or
        constant) in the given ``interval``, False otherwise.

    Examples
    ========

    >>> from sympy import is_decreasing
    >>> from sympy.abc import x, y
    >>> from sympy import S, Interval, oo
    >>> is_decreasing(1/(x**2 - 3*x), Interval.open(S(3)/2, 3))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, S(3)/2))
    False
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, 1.5))
    False
    >>> is_decreasing(-x**2, Interval(-oo, 0))
    False
    >>> is_decreasing(-x**2 + y, Interval(-oo, 0), x)
    False

    """
    return monotonicity_helper(expression, lambda x: x <= 0, interval, symbol)


def is_strictly_decreasing(
    expression: Expr,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Return whether the function is strictly decreasing in the given interval.

    Parameters
    ==========

    expression : Expr
        The target function which is being checked.
    interval : Set, optional
        The range of values in which we are testing (defaults to set of
        all real numbers).
    symbol : Symbol, optional
        The symbol present in expression which gets varied over the given range.

    Returns
    =======

    Boolean
        True if ``expression`` is strictly decreasing in the given ``interval``,
        False otherwise.

    Examples
    ========

    >>> from sympy import is_strictly_decreasing
    >>> from sympy.abc import x, y
    >>> from sympy import S, Interval, oo
    >>> is_strictly_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_strictly_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, S(3)/2))
    False
    >>> is_strictly_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, 1.5))
    False
    >>> is_strictly_decreasing(-x**2, Interval(-oo, 0))
    False
    >>> is_strictly_decreasing(-x**2 + y, Interval(-oo, 0), x)
    False

    """
    return monotonicity_helper(expression, lambda x: x < 0, interval, symbol)


def is_monotonic(
    expression: Expr,
    interval: Set = S.Reals,
    symbol: Symbol | None = None,
) -> bool:
    """
    Return whether the function is monotonic in the given interval.

    Parameters
    ==========

    expression : Expr
        The target function which is being checked.
    interval : Set, optional
        The range of values in which we are testing (defaults to set of
        all real numbers).
    symbol : Symbol, optional
        The symbol present in expression which gets varied over the given range.

    Returns
    =======

    Boolean
        True if ``expression`` is monotonic in the given ``interval``,
        False otherwise.

    Raises
    ======

    NotImplementedError
        Monotonicity check has not been implemented for the queried function.

    Explanation
    ===========

    A real-valued function f is monotonic on an interval I if it is either
    monotonically increasing or monotonically decreasing throughout I:

    - **Monotonically increasing**: For all x1, x2 in I where x1 < x2,
      we have f(x1) <= f(x2)
    - **Monotonically decreasing**: For all x1, x2 in I where x1 < x2,
      we have f(x1) >= f(x2)

    According to Rudin (Theorem 4.30), monotonic functions can have jump
    discontinuities (first kind) but NOT poles (second kind). Poles violate
    the ordering property. Example: tan(x) on [0,5] has poles at pi/2, 3*pi/2
    causing transitions between +-infinity, which breaks monotonicity.

    Examples
    ========

    Basic monotonicity tests:

    >>> from sympy import is_monotonic, tan, Interval, pi, oo
    >>> from sympy.abc import x, y
    >>> from sympy import S
    >>> is_monotonic(x**3 - 3*x**2 + 4*x, S.Reals)
    True
    >>> is_monotonic(-x**2, S.Reals)
    False

    Functions with poles (interior singularities) are not monotonic:

    >>> is_monotonic(tan(x), Interval(0, 5), x)
    False
    >>> is_monotonic(1/x, Interval(-1, 1), x)
    False

    But the same functions are monotonic on intervals that exclude the poles:

    >>> is_monotonic(tan(x), Interval.open(-pi/2, pi/2), x)
    True
    >>> is_monotonic(1/x, Interval.open(0, 1), x)
    True

    Piecewise constant intervals:

    >>> is_monotonic(1/(x**2 - 3*x), Interval.open(S(3)/2, 3), x)
    True
    >>> is_monotonic(1/(x**2 - 3*x), Interval.Lopen(3, oo), x)
    True

    Multivariate expressions:

    >>> is_monotonic(x**2 + y + 1, Interval(1, 2), x)
    True

    See Also
    ========

    is_increasing : Test if function is non-decreasing
    is_decreasing : Test if function is non-increasing
    is_strictly_increasing : Test if function is strictly increasing
    is_strictly_decreasing : Test if function is strictly decreasing
    singularities : Find singularities of a function
    continuous_domain : Find the domain where function is continuous

    References
    ==========

    .. [1] Rudin, W. (1976). "Principles of Mathematical Analysis" (3rd ed.),
           McGraw-Hill. Theorem 4.30: A monotonic function has at most
           countably many discontinuities, and all are of the first kind
           (jump discontinuities).

    .. [2] https://en.wikipedia.org/wiki/Monotonic_function

    """
    from sympy.solvers.solveset import solveset

    expression = sympify(expression)

    free = expression.free_symbols
    if symbol is None and len(free) > 1:
        raise NotImplementedError(
            'is_monotonic has not yet been implemented'
            ' for all multivariate expressions.'
        )

    variable = symbol or (free.pop() if free else Symbol('x'))

    try:
        sings = singularities(expression, variable, interval)
        if interval.is_subset(S.Reals):
            interior_sings = interval.interior.intersection(sings)
            if interior_sings != S.EmptySet:
                return False
    except (NotImplementedError, AttributeError):
        pass

    turning_points = solveset(expression.diff(variable), variable, interval)
    return interval.intersection(turning_points) is S.EmptySet
