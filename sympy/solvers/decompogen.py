from sympy.core import Function, Pow, sympify
from sympy.polys import Poly, decompose


def decompogen(*args):
    """
    Computes General functional decomposition of ``f``.
    Given an expression ``f``, returns a list ``[f_1, f_2, ..., f_n]``,
    where::
              f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n))

    Note: This is a General decomposition function. It also decomposes
    Polynomials. For only Polynomial decomposition see ``decompose`` in polys.

    Examples
    ========

    >>> from sympy.solvers.decompogen import decompogen
    >>> from sympy.abc import x
    >>> from sympy import sqrt, sin, cos
    >>> from sympy.core import Lambda
    >>> decompogen(sin(cos(x)), x)
    [sin(x), cos(x)]
    >>> decompogen(sin(x)**2 + sin(x) + 1, x)
    [x**2 + x + 1, sin(x)]
    >>> decompogen(sqrt(6*x**2 - 5), x)
    [sqrt(x), 6*x**2 - 5]
    >>> decompogen(sin(sqrt(cos(x**2 + 1))), x)
    [sin(x), sqrt(x), cos(x), x**2 + 1]
    >>> decompogen(x**4 + 2*x**3 - x - 1, x)
    [x**2 - x - 1, x**2 + x]
    >>> decompogen(Lambda(x, sin(cos(x))), x)
    [sin(x), cos(x)]
    >>> decompogen(lambda x : sin(cos(x)), x)
    [sin(x), cos(x), x]

    """
    from sympy.core import Lambda
    from sympy.geometry.util import _uniquely_named_symbol
    from sympy.core.function import FunctionClass
    from sympy.utilities.misc import func_name
    from sympy.core.symbol import Symbol, Dummy

    if len(args) is not 2:
        raise ValueError('decompogen expects 2 args (function, symbol), got: %s' % len(args))

    f = args[0]
    symbol = args[1]

    # ===== Lambda Functions ===== #

    if isinstance(f, Lambda):
        f_expr = f.args[1]
        f_x = f.args[0]
        return decompogen(lambda f_x : f_expr, symbol)
    elif (isinstance(f, FunctionClass) or func_name(f) == '<lambda>'):
        var = _uniquely_named_symbol(Symbol('x'), f(Dummy()))
        expr = f(var)
        return decompogen(expr, symbol)

    f = sympify(f)
    result = []

    # ===== Simple Functions ===== #
    if isinstance(f, (Function, Pow)):
        if f.args[0] == symbol:
            return [f]
        result += [f.subs(f.args[0], symbol)] + decompogen(f.args[0], symbol)
        return result

    # ===== Convert to Polynomial ===== #
    fp = Poly(f)
    gens = list(filter(lambda x: symbol in x.free_symbols , fp.gens))

    if len(gens) == 1 and gens[0] != symbol:
        f1 = f.subs(gens[0], symbol)
        f2 = gens[0]
        result += [f1] + decompogen(f2, symbol)
        return result

    # ===== Polynomial decompose() ====== #
    try:
        result += decompose(f)
        return result
    except ValueError:
        return [f]
