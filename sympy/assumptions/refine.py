from sympy.core import S, Add
from sympy.assumptions import Q, ask
from sympy.logic.boolalg import fuzzy_not

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
    if not expr.is_Atom:
        args = [refine(arg, assumptions) for arg in expr.args]
        # TODO: this will probably not work with Integral or Polynomial
        expr = expr.func(*args)
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

    Examples
    ========

    >>> from sympy import Symbol, Q, refine, Abs
    >>> from sympy.assumptions.refine import refine_abs
    >>> from sympy.abc import x
    >>> refine_abs(Abs(x), Q.real(x))
    >>> refine_abs(Abs(x), Q.positive(x))
    x
    >>> refine_abs(Abs(x), Q.negative(x))
    -x

    """
    arg = expr.args[0]
    if ask(Q.real(arg), assumptions) and \
            fuzzy_not(ask(Q.negative(arg), assumptions)):
        # if it's nonnegative
        return arg
    if ask(Q.negative(arg), assumptions):
        return -arg

def refine_Pow(expr, assumptions):
    """
    Handler for instances of Pow.

    >>> from sympy import Symbol, Q
    >>> from sympy.assumptions.refine import refine_Pow
    >>> from sympy.abc import x,y,z
    >>> refine_Pow((-1)**x, Q.real(x))
    >>> refine_Pow((-1)**x, Q.even(x))
    1
    >>> refine_Pow((-1)**x, Q.odd(x))
    -1

    For powers of -1, even parts of the exponent can be simplified:

    >>> refine_Pow((-1)**(x+y), Q.even(x))
    (-1)**y
    >>> refine_Pow((-1)**(x+y+z), Q.odd(x) & Q.odd(z))
    (-1)**y
    >>> refine_Pow((-1)**(x+y+2), Q.odd(x))
    (-1)**(y + 1)
    >>> refine_Pow((-1)**(x+3), True)
    (-1)**(x + 1)

    """
    from sympy.core import Pow, Rational
    from sympy.functions import sign
    if ask(Q.real(expr.base), assumptions):
        if expr.base.is_number:
            if ask(Q.even(expr.exp), assumptions):
                return abs(expr.base) ** expr.exp
            if ask(Q.odd(expr.exp), assumptions):
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

                coeff, terms = expr.exp.as_coeff_add()
                terms = set(terms)
                even_terms = set([])
                odd_terms = set([])
                initial_number_of_terms = len(terms)

                for t in terms:
                    if ask(Q.even(t), assumptions):
                        even_terms.add(t)
                    elif ask(Q.odd(t), assumptions):
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

    >>> from sympy import Symbol, Q, exp, I, pi
    >>> from sympy.assumptions.refine import refine_exp
    >>> from sympy.abc import x
    >>> refine_exp(exp(pi*I*2*x), Q.real(x))
    >>> refine_exp(exp(pi*I*2*x), Q.integer(x))
    1

    """
    arg = expr.args[0]
    if arg.is_Mul:
        coeff = arg.as_coefficient(S.Pi*S.ImaginaryUnit)
        if coeff:
            if ask(Q.integer(2*coeff), assumptions):
                if ask(Q.even(coeff), assumptions):
                    return S.One
                elif ask(Q.odd(coeff), assumptions):
                    return S.NegativeOne
                elif ask(Q.even(coeff + S.Half), assumptions):
                    return -S.ImaginaryUnit
                elif ask(Q.odd(coeff + S.Half), assumptions):
                    return S.ImaginaryUnit

handlers_dict = {
    'Abs'        : refine_abs,
    'Pow'        : refine_Pow,
    'exp'        : refine_exp,
}
