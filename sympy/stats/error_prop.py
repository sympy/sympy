"""Tools for arithmetic error propogation."""
from __future__ import print_function, division
from itertools import repeat

from sympy import Symbol, Add, Mul, simplify, Pow, exp
from sympy.stats.symbolic_probability import RandomSymbol, Variance, Covariance


def variance_prop(expr, consts=()):
    """Symbolically propagates variance (`\sigma^2`) for expressions.
    This is computed as as seen in [1]_.

    Parameters
    ==========
    expr : Expr
        A sympy expression to compute the variance for.
    consts : sequence of Symbols, optional
        Represents symbols that are known constants in the expr,
        and thus have zero variance. All symbols not in consts are
        assumed to be variant.

    Returns
    =======
    var_expr : Expr
        An expression for the total variance of the expr.
        The variance for the original symbols (e.g. x) are represented
        via instance of the Variance symbol (e.g. Variance(x)).

    Examples
    ========

    >>> from sympy import symbols, exp
    >>> from sympy.stats.error_prop import variance_prop
    >>> x, y = symbols('x y')

    >>> variance_prop(x + y)
    Variance(x) + Variance(y)

    >>> variance_prop(x * y)
    Variance(x)*y**2 + Variance(y)*x**2

    >>> variance_prop(exp(2*x))
    4*Variance(x)*exp(4*x)

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    """
    args = expr.args
    if len(args) == 0:
        if expr in consts:
            return 0
        elif isinstance(expr, RandomSymbol):
            return Variance(expr)
        elif isinstance(expr, Symbol):
            return Variance(RandomSymbol(expr))
        else:
            return 0
    var_args = list(map(variance_prop, args, repeat(consts)))
    if isinstance(expr, Add):
        return Add(*var_args)
    elif isinstance(expr, Mul):
        terms = [v/a**2 for a, v in zip(args, var_args)]
        return simplify(expr**2 * Add(*terms))
    elif isinstance(expr, Pow):
        b = args[1]
        v = var_args[0] * (expr * b / args[0])**2
        return simplify(v)
    elif isinstance(expr, exp):
        return simplify(var_args[0] * expr**2)
    else:
        # unknown how to proceed, return variance of whole expr.
        return Variance(expr)


def covariance_prop(expr, consts=()):
    """Symbolically propagates covariance (`\sigma_{AB}`) for expressions.
    This is computed as as seen in [1]_.

    Parameters
    ==========
    expr : Expr
        A sympy expression to compute the covariance for.
    consts : sequence of Symbols, optional
        Represents symbols that are known constants in the expr,
        and thus have zero variance. All symbols not in consts are
        assumed to be variant.

    Returns
    =======
    var_expr : Expr
        An expression for the total variance of the expr.
        The variance for the original symbols (e.g. x) are represented
        via instance of the Variance symbol (e.g. Variance(x)).

    Examples
    ========

    >>> from sympy import symbols, exp
    >>> from sympy.stats.error_prop import covariance_prop
    >>> x, y = symbols('x y')

    >>> covariance_prop(x + y)
    Variance(x) + Variance(y)

    >>> variance_prop(x * y)
    Variance(x)*y**2 + Variance(y)*x**2

    >>> variance_prop(exp(2*x))
    4*Variance(x)*exp(4*x)

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    """
    args = expr.args
    if len(args) == 0:
        if expr in consts:
            return 0
        elif isinstance(expr, RandomSymbol):
            return expr
        elif isinstance(expr, Symbol):
            return RandomSymbol(expr)
        else:
            return 0
    covar_args = list(map(covariance_prop, args, repeat(consts)))
    if isinstance(expr, Add):
        return Add(*var_args)
    elif isinstance(expr, Mul):
        terms = [v/a**2 for a, v in zip(args, var_args)]
        return simplify(expr**2 * Add(*terms))
    elif isinstance(expr, Pow):
        b = args[1]
        v = var_args[0] * (expr * b / args[0])**2
        return simplify(v)
    elif isinstance(expr, exp):
        return simplify(var_args[0] * expr**2)
    else:
        # unknown how to proceed, return variance of whole expr.
        return Variance(expr)
