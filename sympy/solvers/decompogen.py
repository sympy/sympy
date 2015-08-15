from sympy.core import Function, Pow, sympify
from sympy.polys import Poly


def decompogen(f, symbol):
    """
    Computes General functional decomposition of ``f``.
    Given an expression ``f``, returns a list ``[f_1, f_2, ..., f_n]``,
    where::
              f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n))

    Note: This is a General decomposition function. For Polynomial
    decomposition see ``decompose`` in polys.

    Examples
    ========

    >>> from sympy.solvers.decompogen import decompogen
    >>> from sympy.abc import x
    >>> from sympy import sqrt, sin, cos
    >>> decompogen(sin(cos(x)), x)
    [sin(x), cos(x)]
    >>> decompogen(sin(x)**2 + sin(x) + 1, x)
    [x**2 + x + 1, sin(x)]
    >>> decompogen(sqrt(6*x**2 - 5), x)
    [sqrt(x), 6*x**2 - 5]

    """
    f = sympify(f)

    # ===== Simple Functions ===== #
    if isinstance(f, (Function, Pow)):
        return [f.subs(f.args[0], symbol), f.args[0]]

    # ===== Convert to Polynomial ===== #
    fp = Poly(f)
    gens = fp.gens

    if len(gens) == 1:
        f1 = f.subs(gens[0], symbol)
        f2 = gens[0]
        return [f1, f2]

    return [f]
