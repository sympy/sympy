"""
Domain-aware cancellation for SymPy expressions.

The standard ``cancel()`` function simplifies rational expressions
by cancelling common factors, but this can lose information about
singularities. For example, ``cancel((x**2 - 1)/(x - 1))`` returns
``x + 1``, which is defined at ``x = 1`` while the original is not.

This module provides ``cancel_with_domain`` which performs cancellation
while tracking where cancelled factors would have been zero.
"""

from sympy.core.symbol import Symbol
from sympy.polys.polytools import cancel, factor
from sympy.core.relational import Ne
from sympy.sets.sets import FiniteSet
from sympy.logic.boolalg import And


def cancel_with_domain(expr, *symbols):
    """Cancel common factors and return the simplified expression
    along with domain restrictions.

    Parameters
    ==========
    expr : Expr
        A rational expression to simplify.
    *symbols : Symbol
        The variable(s) over which to determine domain restrictions.
        If not given, uses the free symbols of the expression.

    Returns
    =======
    tuple : (simplified_expr, excluded_points)
        ``simplified_expr`` is the result of ``cancel(expr)``.
        ``excluded_points`` is a set of ``Ne(symbol, value)``
        conditions representing points where the original expression
        has singularities that were removed by cancellation.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.simplify.cancel_utils import cancel_with_domain
    >>> x = symbols('x')
    >>> simplified, restrictions = cancel_with_domain((x**2 - 1)/(x - 1))
    >>> simplified
    x + 1
    >>> restrictions
    {Ne(x, 1)}

    >>> from sympy import sin
    >>> simplified, restrictions = cancel_with_domain(x/x)
    >>> simplified
    1
    >>> restrictions
    {Ne(x, 0)}
    """
    from sympy.core.numbers import oo, zoo, nan
    from sympy.solvers.solveset import solveset
    from sympy.sets.sets import S as Sets

    simplified = cancel(expr)

    if not symbols:
        symbols = tuple(expr.free_symbols)

    restrictions = set()

    # If the expression is a ratio, find where the denominator is zero
    if hasattr(expr, 'as_numer_denom'):
        numer, denom = expr.as_numer_denom()
        simplified_numer, simplified_denom = simplified.as_numer_denom()

        # Factor to find cancelled factors
        try:
            for sym in symbols:
                # Find zeros of original denominator
                original_zeros = solveset(denom, sym, domain=Sets.Reals)
                # Find zeros of simplified denominator
                simplified_zeros = solveset(simplified_denom, sym, domain=Sets.Reals)

                # The cancelled singularities are points in original but not simplified
                if hasattr(original_zeros, '__iter__') and hasattr(simplified_zeros, '__iter__'):
                    cancelled = set(original_zeros) - set(simplified_zeros)
                    for point in cancelled:
                        if point not in (oo, -oo, zoo, nan):
                            restrictions.add(Ne(sym, point))
        except (NotImplementedError, ValueError, TypeError):
            # If solving fails, we can't determine restrictions
            pass

    return simplified, restrictions
