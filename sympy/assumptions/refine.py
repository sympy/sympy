from __future__ import print_function, division

from sympy.core import S, Add, Expr, Basic
from sympy.assumptions import Q, ask, assuming
from sympy.core.logic import fuzzy_not


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
    name = expr.__class__.__name__
    handler = handlers_dict.get(name, None)
    if handler is None:
        return expr
    with assuming(assumptions):
        new_expr = handler(expr)
    if (new_expr is None) or (expr == new_expr):
        return expr
    if not isinstance(new_expr, Expr):
        return new_expr
    return refine(new_expr, assumptions)


def refine_abs(expr):
    """
    Handler for the absolute value.

    Examples
    ========

    >>> from sympy import Symbol, Q, refine, Abs
    >>> from sympy.assumptions.refine import refine
    >>> from sympy.abc import x
    >>> refine(Abs(x), Q.real(x))
    >>> refine(Abs(x), Q.positive(x))
    x
    >>> refine(Abs(x), Q.negative(x))
    -x

    """
    arg = expr.args[0]
    if ask(Q.real(arg)) and \
            fuzzy_not(ask(Q.negative(arg))):
        # if it's nonnegative
        return arg
    if ask(Q.negative(arg)):
        return -arg


def refine_Pow(expr):
    """
    Handler for instances of Pow.

    >>> from sympy import Symbol, Q
    >>> from sympy.assumptions.refine import refine
    >>> from sympy.abc import x,y,z
    >>> refine((-1)**x, Q.real(x))
    >>> refine((-1)**x, Q.even(x))
    1
    >>> refine((-1)**x, Q.odd(x))
    -1

    For powers of -1, even parts of the exponent can be simplified:

    >>> refine((-1)**(x+y), Q.even(x))
    (-1)**y
    >>> refine((-1)**(x+y+z), Q.odd(x) & Q.odd(z))
    (-1)**y
    >>> refine((-1)**(x+y+2), Q.odd(x))
    (-1)**(y + 1)
    >>> refine((-1)**(x+3), True)
    (-1)**(x + 1)

    """
    from sympy.core import Pow, Rational
    from sympy.functions.elementary.complexes import Abs
    from sympy.functions import sign
    if isinstance(expr.base, Abs):
        if ask(Q.real(expr.base.args[0])) and \
                ask(Q.even(expr.exp)):
            return expr.base.args[0] ** expr.exp
    if ask(Q.real(expr.base)):
        if expr.base.is_number:
            if ask(Q.even(expr.exp)):
                return abs(expr.base) ** expr.exp
            if ask(Q.odd(expr.exp)):
                return sign(expr.base) * abs(expr.base) ** expr.exp
        if isinstance(expr.exp, Rational):
            if type(expr.base) is Pow:
                return abs(expr.base.base) ** (expr.base.exp * expr.exp)

        if expr.base is S.NegativeOne:
            if expr.exp.is_Add:

                old = expr

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
                    if ask(Q.even(t)):
                        even_terms.add(t)
                    elif ask(Q.odd(t)):
                        odd_terms.add(t)

                terms -= even_terms
                if len(odd_terms) % 2:
                    terms -= odd_terms
                    new_coeff = (coeff + S.One) % 2
                else:
                    terms -= odd_terms
                    new_coeff = coeff % 2

                if new_coeff != coeff or len(terms) < initial_number_of_terms:
                    terms.add(new_coeff)
                    expr = expr.base**(Add(*terms))

                # Handle (-1)**((-1)**n/2 + m/2)
                e2 = 2*expr.exp
                if ask(Q.even(e2)):
                    if e2.could_extract_minus_sign():
                        e2 *= expr.base
                if e2.is_Add:
                    i, p = e2.as_two_terms()
                    if p.is_Pow and p.base is S.NegativeOne:
                        if ask(Q.integer(p.exp)):
                            i = (i + 1)/2
                            if ask(Q.even(i)):
                                return expr.base**p.exp
                            elif ask(Q.odd(i)):
                                return expr.base**(p.exp + 1)
                            else:
                                return expr.base**(p.exp + i)

                if old != expr:
                    return expr


def refine_exp(expr):
    """
    Handler for exponential function.

    >>> from sympy import Symbol, Q, exp, I, pi
    >>> from sympy.assumptions.refine import refine
    >>> from sympy.abc import x
    >>> refine(exp(pi*I*2*x), Q.real(x))
    >>> refine(exp(pi*I*2*x), Q.integer(x))
    1

    """
    arg = expr.args[0]
    if arg.is_Mul:
        coeff = arg.as_coefficient(S.Pi*S.ImaginaryUnit)
        if coeff:
            if ask(Q.integer(2*coeff)):
                if ask(Q.even(coeff)):
                    return S.One
                elif ask(Q.odd(coeff)):
                    return S.NegativeOne
                elif ask(Q.even(coeff + S.Half)):
                    return -S.ImaginaryUnit
                elif ask(Q.odd(coeff + S.Half)):
                    return S.ImaginaryUnit

def refine_Relational(expr):
    """
    Handler for Relational

    >>> from sympy.assumptions.refine import refine
    >>> from sympy.assumptions.ask import Q
    >>> from sympy.abc import x
    >>> refine(x<0, ~Q.is_true(x<0))
    False
    """
    return ask(Q.is_true(expr))


handlers_dict = {
    'Abs': refine_abs,
    'Pow': refine_Pow,
    'exp': refine_exp,
    'Equality' : refine_Relational,
    'Unequality' : refine_Relational,
    'GreaterThan' : refine_Relational,
    'LessThan' : refine_Relational,
    'StrictGreaterThan' : refine_Relational,
    'StrictLessThan' : refine_Relational
}
