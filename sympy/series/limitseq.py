"""Limits of sequences"""

from __future__ import print_function, division

from sympy.core.sympify import sympify
from sympy.core.singleton import S


def difference_delta(expr, n=None, step=1):
    """Difference Operator.

    Discrete analogous to differential operator.

    Examples
    ========

    >>> from sympy import difference_delta as dd
    >>> from sympy.abc import n
    >>> dd(n*(n + 1), n)
    2*n + 2
    >>> dd(n*(n + 1), n, 2)
    4*n + 6

    References
    ==========

    .. [1] https://reference.wolfram.com/language/ref/DifferenceDelta.html
    """
    expr = sympify(expr)

    if n is None:
        f = expr.free_symbols
        if len(f) == 1:
            n = f.pop()
        elif len(f) == 0:
            return S.Zero
        else:
            raise ValueError("Since there is more than one variable in the"
                             " expression, a variable must be supplied to"
                             " take the difference of %s" % expr)
    step = sympify(step)
    if step.is_number is False:
        raise ValueError("Step should be a number.")
    elif step in [S.Infinity, -S.Infinity]:
        raise ValueError("Step should be bounded.")

    if hasattr(expr, '_eval_difference_delta'):
        return expr._eval_difference_delta(n, step)

    return expr.subs(n, n + step) - expr
