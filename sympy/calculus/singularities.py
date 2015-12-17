from sympy.core.sympify import sympify
from sympy.solvers.solveset import solveset
from sympy.simplify import simplify
from sympy import S


def singularities(expr, sym):
    """
    Finds singularities for a function.
    Currently supported functions are:
    - univariate real rational functions

    Examples
    ========

    >>> from sympy.calculus.singularities import singularities
    >>> from sympy import Symbol
    >>> x = Symbol('x', real=True)
    >>> singularities(x**2 + x + 1, x)
    ()
    >>> singularities(1/(x + 1), x)
    (-1,)

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Mathematical_singularity

    """
    if not expr.is_rational_function(sym):
        raise NotImplementedError("Algorithms finding singularities for"
                                  " non rational functions are not yet"
                                  " implemented")
    else:
        return solveset(simplify(1/expr), sym)

###########################################################################
###################### DIFFERENTIAL CALCULUS METHODS ######################
###########################################################################


def is_increasing(f, interval=S.Reals, symbol=None):
    """
    Returns if a function is increasing or not, in the given
    ``Interval``.

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
    f = sympify(f)
    free_sym = f.free_symbols

    if symbol is None:
        if len(free_sym) > 1:
            raise NotImplementedError('is_increasing has not yet been implemented '
                                        'for all multivariate expressions')
        if len(free_sym) == 0:
            return True
        symbol = free_sym.pop()

    df = f.diff(symbol)
    df_nonneg_interval = solveset(df >= 0, symbol, domain=S.Reals)
    return interval.is_subset(df_nonneg_interval)


def is_strictly_increasing(f, interval=S.Reals, symbol=None):
    """
    Returns if a function is strictly increasing or not, in the given
    ``Interval``.

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
    f = sympify(f)
    free_sym = f.free_symbols

    if symbol is None:
        if len(free_sym) > 1:
            raise NotImplementedError('is_strictly_increasing has not yet been implemented '
                                        'for all multivariate expressions')
        elif len(free_sym) == 0:
            return False
        symbol = free_sym.pop()

    df = f.diff(symbol)
    df_pos_interval = solveset(df > 0, symbol, domain=S.Reals)
    return interval.is_subset(df_pos_interval)


def is_decreasing(f, interval=S.Reals, symbol=None):
    """
    Returns if a function is decreasing or not, in the given
    ``Interval``.

    Examples
    ========

    >>> from sympy import is_decreasing
    >>> from sympy.abc import x, y
    >>> from sympy import S, Interval, oo
    >>> is_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, S(3)/2))
    False
    >>> is_decreasing(-x**2, Interval(-oo, 0))
    False
    >>> is_decreasing(-x**2 + y, Interval(-oo, 0), x)
    False

    """
    f = sympify(f)
    free_sym = f.free_symbols

    if symbol is None:
        if len(free_sym) > 1:
            raise NotImplementedError('is_decreasing has not yet been implemented '
                                        'for all multivariate expressions')
        elif len(free_sym) == 0:
            return True
        symbol = free_sym.pop()

    df = f.diff(symbol)
    df_nonpos_interval = solveset(df <= 0, symbol, domain=S.Reals)
    return interval.is_subset(df_nonpos_interval)


def is_strictly_decreasing(f, interval=S.Reals, symbol=None):
    """
    Returns if a function is strictly decreasing or not, in the given
    ``Interval``.

    Examples
    ========

    >>> from sympy import is_strictly_decreasing
    >>> from sympy.abc import x, y
    >>> from sympy import S, Interval, oo
    >>> is_strictly_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))
    True
    >>> is_strictly_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_strictly_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, S(3)/2))
    False
    >>> is_strictly_decreasing(-x**2, Interval(-oo, 0))
    False
    >>> is_strictly_decreasing(-x**2 + y, Interval(-oo, 0), x)
    False

    """
    f = sympify(f)
    free_sym = f.free_symbols

    if symbol is None:
        if len(free_sym) > 1:
            raise NotImplementedError('is_strictly_decreasing has not yet been implemented '
                                        'for all multivariate expressions')
        elif len(free_sym) == 0:
            return False
        symbol = free_sym.pop()

    df = f.diff(symbol)
    df_neg_interval = solveset(df < 0, symbol, domain=S.Reals)
    return interval.is_subset(df_neg_interval)


def is_monotonic(f, interval=S.Reals, symbol=None):
    """
    Returns if a function is monotonic or not, in the given
    ``Interval``.

    Examples
    ========

    >>> from sympy import is_monotonic
    >>> from sympy.abc import x, y
    >>> from sympy import S, Interval, oo
    >>> is_monotonic(1/(x**2 - 3*x), Interval.open(1.5, 3))
    True
    >>> is_monotonic(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_monotonic(x**3 - 3*x**2 + 4*x, S.Reals)
    True
    >>> is_monotonic(-x**2, S.Reals)
    False
    >>> is_monotonic(x**2 + y + 1, Interval(1, 2), x)
    True

    """
    from sympy.core.logic import fuzzy_or
    f = sympify(f)
    free_sym = f.free_symbols

    if symbol is None and len(free_sym) > 1:
        raise NotImplementedError('is_monotonic has not yet been '
                                'for all multivariate expressions')

    inc = is_increasing(f, interval, symbol)
    dec = is_decreasing(f, interval, symbol)

    return fuzzy_or([inc, dec])
