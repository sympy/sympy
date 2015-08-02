from sympy.solvers import solve
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
        return tuple(sorted(solveset(simplify(1/expr), sym)))

###########################################################################
###################### DIFFERENTIAL CALCULUS METHODS ######################
###########################################################################


def is_increasing(f, given_interval):
    """
    Returns if a function is increasing or not, in the given
    `Interval`.

    Examples
    ========

    >>> from sympy.calculus.singularities import is_increasing
    >>> from sympy.abc import x
    >>> from sympy import S, Interval, oo
    >>> is_increasing(x**3 - 3*x**2 + 4*x, S.Reals)
    True
    >>> is_increasing(-x**2, Interval(-oo, 0))
    True
    >>> is_increasing(-x**2, Interval(0, oo))
    False
    >>> is_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval(-2, 3))
    False

    """
    if len(f.free_symbols) > 1:
        raise NotImplementedError
    symbol = f.free_symbols.pop()
    df = f.diff(symbol)
    df_pos_interval = solveset(df >= 0, symbol, domain=S.Reals)
    return given_interval.is_subset(df_pos_interval)


def is_strictly_increasing(f, given_interval):
    """
    Returns if a function is strictly increasing or not, in the given
    `Interval`.

    Examples
    ========

    >>> from sympy.calculus.singularities import is_strictly_increasing
    >>> from sympy.abc import x
    >>> from sympy import S, Interval, oo
    >>> is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.Ropen(-oo, -2))
    True
    >>> is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.Lopen(3, oo))
    True
    >>> is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.open(-2, 3))
    False
    >>> is_strictly_increasing(-x**2, Interval(0, oo))
    False

    """
    if len(f.free_symbols) > 1:
        raise NotImplementedError
    symbol = f.free_symbols.pop()
    df = f.diff(symbol)
    df_pos_interval = solveset(df > 0, symbol, domain=S.Reals)
    return given_interval.is_subset(df_pos_interval)


def is_decreasing(f, given_interval):
    """
    Returns if a function is decreasing or not, in the given
    `Interval`.

    Examples
    ========

    >>> from sympy.calculus.singularities import is_decreasing
    >>> from sympy.abc import x
    >>> from sympy import S, Interval, oo
    >>> is_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, S(3)/2))
    False
    >>> is_decreasing(-x**2, Interval(-oo, 0))
    False

    """
    if len(f.free_symbols) > 1:
        raise NotImplementedError
    symbol = f.free_symbols.pop()
    df = f.diff(symbol)
    df_nonpos_interval = solveset(df <= 0, symbol, domain=S.Reals)
    return given_interval.is_subset(df_nonpos_interval)


def is_strictly_decreasing(f, given_interval):
    """
    Returns if a function is decreasing or not, in the given
    `Interval`.

    Examples
    ========

    >>> from sympy.calculus.singularities import is_decreasing
    >>> from sympy.abc import x
    >>> from sympy import S, Interval, oo
    >>> is_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    True
    >>> is_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, S(3)/2))
    False
    >>> is_decreasing(-x**2, Interval(-oo, 0))
    False

    """
    if len(f.free_symbols) > 1:
        raise NotImplementedError
    symbol = f.free_symbols.pop()
    df = f.diff(symbol)
    df_neg_interval = solveset(df < 0, symbol, domain=S.Reals)
    return given_interval.is_subset(df_neg_interval)
