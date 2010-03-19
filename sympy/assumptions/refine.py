from sympy.core import S, Add
from sympy.utilities.source import get_class
from sympy.assumptions import Q, ask
from sympy.logic.boolalg import fuzzy_not

def refine(expr, assumptions=True):
    """
    Simplify an expression using assumptions.

    Gives the form of expr that would be obtained if symbols
    in it were replaced by explicit numerical expressions satisfying
    the assumptions.

    Examples::

        >>> from sympy import refine, sqrt, Assume, Q
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

    Examples::

    >>> from sympy import Symbol, Assume, Q, refine
    >>> from sympy.assumptions.refine import refine_abs
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
    >>> from sympy.assumptions.refine import refine_Pow
    >>> from sympy.abc import x,y,z
    >>> refine_Pow((-1)**x, Assume(x, Q.real))
    >>> refine_Pow((-1)**x, Assume(x, Q.even))
    1
    >>> refine_Pow((-1)**x, Assume(x, Q.odd))
    -1

    For powers of -1, even parts of the exponent can be simplified:

    >>> refine_Pow((-1)**(x+y), Assume(x, Q.even))
    (-1)**y
    >>> refine_Pow((-1)**(x+y+z), Assume(x, Q.odd) & Assume(z, Q.odd))
    (-1)**y
    >>> refine_Pow((-1)**(x+y+2), Assume(x, Q.odd))
    (-1)**(1 + y)
    >>> refine_Pow((-1)**(x+3), True)
    (-1)**(1 + x)


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

        if expr.base is S.NegativeOne:
            if expr.exp.is_Add:

                # For powers of (-1) we can remove
                #  - even terms
                #  - pairs of odd terms
                #  - a single odd term + 1
                #  - A numerical constant N can be replaced with mod(N,2)

                coeff, terms = expr.exp.as_coeff_factors()
                terms = set(terms)
                even_terms = set([])
                odd_terms = set([])
                initial_number_of_terms = len(terms)

                for t in terms:
                    if ask(t, Q.even, assumptions):
                        even_terms.add(t)
                    elif ask(t, Q.odd, assumptions):
                        odd_terms.add(t)

                terms -= even_terms
                if len(odd_terms)%2:
                    terms -= odd_terms
                    new_coeff = (coeff + S.One) % 2
                else:
                    terms -= odd_terms
                    new_coeff = coeff % 2

                if new_coeff != coeff or len(terms) < initial_number_of_terms:
                    terms.add(new_coeff)
                    return expr.base**(Add(*terms))


def refine_exp(expr, assumptions):
    """
    Handler for exponential function.

    >>> from sympy import Symbol, Assume, Q, exp, I, pi
    >>> from sympy.assumptions.refine import refine_exp
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
