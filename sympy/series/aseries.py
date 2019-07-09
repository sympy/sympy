from __future__ import print_function, division

from sympy.core.sympify import sympify


def aseries(expr, x=None, n=6, bound=0, hir=False):
    """Asymptotic Series expansion of expr.
    This is equivalent to ``expr.series(x, oo, n)``.

    Parameters
    ==========

    expr : Expression
           The expression whose series is to be expanded.

    x : Symbol
        It is the variable of the expression to be calculated.

    n : Value
        The number of terms upto which the series is to be expanded.

    hir : Boolean
          Set this parameter to be True to produce hierarchical series.
          It stops the recursion at an early level and may provide nicer
          and more useful results.

    bound : Value, Integer
            Use the ``bound`` parameter to give limit on rewriting
            coefficients in its normalised form.

    Examples
    ========

    >>> from sympy import sin, exp
    >>> from sympy.abc import x, y

    >>> e = sin(1/x + exp(-x)) - sin(1/x)

    >>> e.aseries(x)
    (1/(24*x**4) - 1/(2*x**2) + 1 + O(x**(-6), (x, oo)))*exp(-x)

    >>> e.aseries(x, n=3, hir=True)
    -exp(-2*x)*sin(1/x)/2 + exp(-x)*cos(1/x) + O(exp(-3*x), (x, oo))

    >>> e = exp(exp(x)/(1 - 1/x))

    >>> e.aseries(x, bound=3)
    exp(exp(x)/x**2)*exp(exp(x)/x)*exp(-exp(x) + exp(x)/(1 - 1/x) - exp(x)/x - exp(x)/x**2)*exp(exp(x))

    >>> e.aseries(x)
    exp(exp(x)/(1 - 1/x))

    Returns
    =======

    Expr
        Asymptotic series expansion of the expression.

    Notes
    =====

    This algorithm is directly induced from the limit computational algorithm provided by Gruntz.
    It majorly uses the mrv and rewrite sub-routines. The overall idea of this algorithm is first
    to look for the most rapidly varying subexpression w of a given expression f and then expands f
    in a series in w. Then same thing is recursively done on the leading coefficient
    till we get constant coefficients.

    If the most rapidly varying subexpression of a given expression f is f itself,
    the algorithm tries to find a normalised representation of the mrv set and rewrites f
    using this normalised representation.

    If the expansion contains an order term, it will be either ``O(x**(-n))`` or ``O(w**(-n))``
    where ``w`` belongs to the most rapidly varying expression of ``self``.

    References
    ==========

    .. [1] A New Algorithm for Computing Asymptotic Series - Dominik Gruntz
    .. [2] Gruntz thesis - p90
    .. [3] http://en.wikipedia.org/wiki/Asymptotic_expansion

    See Also
    ========

    See the docstring of Expr.aseries() for complete details of this wrapper.
    """
    expr = sympify(expr)
    return expr.aseries(x, n, bound, hir)