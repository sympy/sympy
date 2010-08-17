"""Tools for simplification of rational expressions. """

from sympy.core import Basic, Add, sympify
from sympy.core.exprtools import gcd_terms

def together(expr, deep=False, symbolic=False):
    """
    Denest and combine rational expressions using symbolic methods.

    This function takes an expression or a container of expressions
    and puts it (them) together by denesting and combining rational
    subexpressions. No heroic measures are taken to minimize degree
    of the resulting numerator and denominator. To obtain completely
    reduced expression use :func:`cancel`. However, :func:`together`
    can preserve as much as possible of the structure of the input
    expression in the output (no expansion is performed).

    A wide variety of objects can be put together including lists,
    tuples, sets, relational objects, integrals and others. It is
    also possible to transform interior of function applications,
    by setting ``deep`` flag to ``True``.

    By definition, :func:`together` is a complement to :func:`apart`,
    so ``apart(together(expr))`` should return expr unchanged. Note
    however, that :func:`together` uses only symbolic methods, so
    it might be necessary to use :func:`cancel` to perform algebraic
    simplification and minimise degree of the numerator and denominator.

    Example
    =======

    >>> from sympy import together, exp
    >>> from sympy.abc import x, y, z

    >>> together(1/x + 1/y)
    (x + y)/(x*y)
    >>> together(1/x + 1/y + 1/z)
    (x*y + x*z + y*z)/(x*y*z)

    >>> together(1/(x*y) + 1/y**2)
    (x + y)/(x*y**2)

    >>> together(1/(1 + 1/x) + 1/(1 + 1/y))
    (x*(1 + y) + y*(1 + x))/((1 + x)*(1 + y))

    >>> together(exp(1/x + 1/y))
    exp(1/x + 1/y)
    >>> together(exp(1/x + 1/y), deep=True)
    exp((x + y)/(x*y))

    >>> together(1/exp(x) + 1/(x*exp(x)))
    (1 + x)*exp(-x)/x

    >>> together(1/exp(2*x) + 1/(x*exp(3*x)))
    (1 + x*exp(x))*exp(-3*x)/x

    """
    def _together(expr):
        if isinstance(expr, Basic):
            if expr.is_commutative is False:
                numer, denom = expr.as_numer_denom()
                return numer/denom
            elif expr.is_Atom or (expr.is_Function and not deep):
                return expr
            elif expr.is_Add:
                return gcd_terms(map(_together, Add.make_args(expr)))
            elif expr.is_Pow:
                base = _together(expr.base)

                if deep:
                    exp = _together(expr.exp)
                else:
                    exp = expr.exp

                return expr.__class__(base, exp)
            else:
                return expr.__class__(*[ _together(arg) for arg in expr.args ])
        elif hasattr(expr, '__iter__'):
            return expr.__class__([ _together(ex) for ex in expr ])

        return expr

    if not symbolic:
        return _together(sympify(expr))
    else:
        from sympy.simplify.simplify import powsimp, separate
        return powsimp(_together(separate(expr)), deep=True, combine='exp')
