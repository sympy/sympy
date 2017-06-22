"""Tools for manipulation of rational expressions. """

from __future__ import print_function, division

from sympy.core import Basic, Add, sympify
from sympy.core.compatibility import iterable
from sympy.core.exprtools import gcd_terms, factor_terms
from sympy.core.function import expand_mul
from sympy.core.mul import Mul
from sympy.utilities import public

@public
def together(expr, deep=False, _expand=None):
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

    Examples
    ========

    >>> from sympy import together, exp
    >>> from sympy.abc import x, y, z

    >>> together(1/x + 1/y)
    (x + y)/(x*y)
    >>> together(1/x + 1/y + 1/z)
    (x*y + x*z + y*z)/(x*y*z)

    >>> together(1/(x*y) + 1/y**2)
    (x + y)/(x*y**2)

    >>> together(1/(1 + 1/x) + 1/(1 + 1/y))
    (x*(y + 1) + y*(x + 1))/((x + 1)*(y + 1))

    >>> together(exp(1/x + 1/y))
    exp(1/y + 1/x)
    >>> together(exp(1/x + 1/y), deep=True)
    exp((x + y)/(x*y))

    >>> together(1/exp(x) + 1/(x*exp(x)))
    (x + 1)*exp(-x)/x

    >>> together(1/exp(2*x) + 1/(x*exp(3*x)))
    (x*exp(x) + 1)*exp(-3*x)/x

    """
    from sympy.simplify.radsimp import fraction
    def _together(expr):
        if isinstance(expr, Basic):
            if expr.is_Atom or (expr.is_Function and not deep):
                return expr
            elif expr.is_Add:
                rv = gcd_terms(list(map(_together, expr.args)))
            elif expr.is_Mul:
                rv = expr.__class__(*[ _together(arg) for arg in expr.args ])
            elif expr.is_Pow:
                base = _together(expr.base)

                if deep:
                    exp = _together(expr.exp)
                else:
                    exp = expr.exp

                return expr.__class__(base, exp)
            else:
                return expr.__class__(*[ _together(arg) for arg in expr.args ])
            if _expand:
                # handle Mul and Add
                old_n, d = fraction(rv, exact=True)
                n = factor_terms(expand_mul(old_n))
                if n.count_ops() >= old_n.count_ops():
                    return rv
                if n.is_Add and d.is_Number:
                    return Mul(1/d, n, evaluate=False)
                elif n != old_n:
                    return n/d
            return rv
        elif iterable(expr):
            return expr.__class__([ _together(ex) for ex in expr ])

        return expr

    return _together(sympify(expr))
