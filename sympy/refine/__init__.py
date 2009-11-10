from sympy.core import S, Symbol, sympify
from sympy.utilities.source import get_class
from sympy.queries import Q, ask
from sympy.logic.boolalg import fuzzy_not

def refine(expr, assumptions=True):
    """
    Simplify an expression using assumptions

    Gives the form of expr that would be obtained if symbols
    in it were replaced by explicit numerical expressions satisfying
    the assumptions.

    Examples:

        >>> from sympy import *
        >>> from sympy.abc import x
        >>> refine(sqrt(x**2), Assume(x, Q.real))
        abs(x)
        >>> refine(sqrt(x**2), Assume(x, Q.positive))
        x

    """
    if not expr.is_Atom:
        args = [refine(arg, assumptions) for arg in expr.args]
        # TODO: this will probably not work with Integral or Polynomial
        expr = type(expr)(*args)
    name = expr.__class__.__name__
    handler = handlers_dict.get(name, None)
    if handler is None: return expr
    new_expr = handler(expr, assumptions)
    if (new_expr is None) or (expr == new_expr):
        return expr
    return refine(new_expr, assumptions)

def refine_abs(expr, assumptions):
    """
    Handler for the absolute value.

    Examples:

    >>> from sympy import Symbol, Assume, Q
    >>> from sympy.refine import refine_abs
    >>> from sympy.abc import x
    >>> refine_abs(abs(x), Assume(x, Q.real))
    >>> refine_abs(abs(x), Assume(x, Q.positive))
    x
    >>> refine_abs(abs(x), Assume(x, Q.negative))
    -x

    """
    arg = expr.args[0]
    if ask(arg, Q.real, assumptions) and \
            fuzzy_not(ask(arg, Q.negative, assumptions)):
        # if it's nonnegative
        return arg
    if ask(arg, Q.negative, assumptions):
        return -arg

def refine_Pow(expr, assumptions):
    """
    Handler for instances of Pow.

    >>> from sympy import Symbol, Assume, Q
    >>> from sympy.refine import refine_Pow
    >>> from sympy.abc import x
    >>> refine_Pow((-1)**x, Assume(x, Q.real))
    >>> refine_Pow((-1)**x, Assume(x, Q.even))
    1
    >>> refine_Pow((-1)**x, Assume(x, Q.odd))
    -1

    """
    from sympy.core import Pow, Rational
    from sympy.functions import sign
    if ask(expr.base, Q.real, assumptions):
        if expr.base.is_number:
            if ask(expr.exp, Q.even, assumptions):
                return abs(expr.base) ** expr.exp
            if ask(expr.exp, Q.odd, assumptions):
                return sign(expr.base) * abs(expr.base) ** expr.exp
        if isinstance(expr.exp, Rational):
            if type(expr.base) is Pow:
                return abs(expr.base.base) ** (expr.base.exp * expr.exp)

def refine_exp(expr, assumptions):
    """
    Handler for exponential function.

    >>> from sympy import Symbol, Assume, Q, exp, I, pi
    >>> from sympy.refine import refine_exp
    >>> from sympy.abc import x
    >>> refine_exp(exp(pi*I*2*x), Assume(x, Q.real))
    >>> refine_exp(exp(pi*I*2*x), Assume(x, Q.integer))
    1

    """
    arg = expr.args[0]
    if arg.is_Mul:
        coeff = arg.as_coefficient(S.Pi*S.ImaginaryUnit)
        if coeff:
            if ask(2*coeff, Q.integer, assumptions):
                if ask(coeff, Q.even, assumptions):
                    return S.One
                elif ask(coeff, Q.odd, assumptions):
                    return S.NegativeOne
                elif ask(coeff + S.Half, Q.even, assumptions):
                    return -S.ImaginaryUnit
                elif ask(coeff + S.Half, Q.odd, assumptions):
                    return S.ImaginaryUnit

handlers_dict = {
    'abs'        : refine_abs,
    'Pow'        : refine_Pow,
    'exp'        : refine_exp,
}
