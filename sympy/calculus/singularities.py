from sympy.solvers import solve
from sympy.simplify import simplify


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
        raise NotImplementedError("Sorry, Algorithms finding singularities for"
                                  " non rational funcitons are not yet"
                                  " implemented")
    else:
        return tuple(sorted(solve(simplify(1/expr), sym)))
