from sympy.polys.polytools import together
from sympy.solvers.solvers import solve, denoms
from sympy.core.sympify import sympify
from sympy.core.symbol import Dummy
from sympy.functions.elementary.exponential import exp


def singularities(expr, sym=None, strict=False, poles=True):
    """Return the real-valued singularities of an expression that
    would set any denominator in the expression to zero.

    Examples
    ========

    >>> from sympy.calculus.singularities import singularities
    >>> from sympy import Symbol, log
    >>> from sympy.abc import x
    >>> singularities(x**2 + x + 1, x)
    ()
    >>> singularities(1/(1/x - 1))
    (1,)
    >>> singularities(1/(1/x - 1), strict=True)
    (0, 1)

    If `poles` is True then points at which functions are undefined will
    be returned. (This is currently not implemented and an error will be
    raised unless the expression is rational.)

    >>> singularities(log(1 - x)/(x - 2), poles=False)
    (2,)
    >>> singularities(log(1 - x)/(x - 2))
    Traceback (most recent call last):
    ...
    NotImplementedError: finding poles for non-rational functions

    Parameters
    ==========

    expr : expression
    sym : the independent symbol (automatically detected for univariate expr)
    strict : if True (default=False) removable singularities are returned, too

    Notes
    =====

    The purpose of this function is to identify any value that would create
    a "division by zero" error. This function does not try to identify the
    domain of the expression, e.g. no singularities for log(x) are identified.

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Mathematical_singularity

    """
    expr = sympify(expr)
    if sym is None:
        free = expr.free_symbols
        if free:
            if len(free) > 1:
                raise ValueError('must supply independent symbol')
            sym = free.pop()

    if not expr.has(sym):
        return tuple()

    if not sym.is_real:
        d = Dummy(sym.name, real=True)
        expr = expr.xreplace({sym: d})
        sym = d

    if not strict:
        expr = together(expr, deep=True)
    rv = set()
    for d in denoms(expr, [sym]):
        rv.update(solve(d, sym))
    # XXX add in poles
    if poles and not expr.is_rational_function(sym):
        raise NotImplementedError('finding poles for non-rational functions')
    return tuple(sorted(rv))
