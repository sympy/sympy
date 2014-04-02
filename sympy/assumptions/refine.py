from __future__ import print_function, division

from sympy.core import S, Add, Expr, Basic
from sympy.assumptions import Q, ask, assuming
from sympy.core.logic import fuzzy_not
from sympy.dispatch import dispatch


def refine(expr, assumptions=True):
    """
    Simplify an expression using assumptions.

    Gives the form of expr that would be obtained if symbols
    in it were replaced by explicit numerical expressions satisfying
    the assumptions.

    Examples
    ========

        >>> from sympy import refine, sqrt, Q
        >>> from sympy.abc import x
        >>> refine(sqrt(x**2), Q.real(x))
        Abs(x)
        >>> refine(sqrt(x**2), Q.positive(x))
        x

    """
    if not isinstance(expr, Basic):
        return expr
    if not expr.is_Atom:
        args = [refine(arg, assumptions) for arg in expr.args]
        # TODO: this will probably not work with Integral or Polynomial
        expr = expr.func(*args)
    with assuming(assumptions):
        new_expr = _refine(expr)
    if (new_expr is None) or (expr == new_expr):
        return expr
    if not isinstance(new_expr, Expr):
        return new_expr
    return refine(new_expr, assumptions)

@dispatch(Basic)
def _refine(expr):
    return expr
