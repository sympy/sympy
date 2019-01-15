from __future__ import print_function, division

from sympy.core.sympify import sympify

from sympy.series.limits import limit
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.function import diff

def series(expr, x=None, x0=0, n=6, dir="+"):
    """Series expansion of expr around point `x = x0`.

    See the docstring of Expr.series() for complete details of this wrapper.
    """
    expr = sympify(expr)
    return expr.series(x, x0, n, dir)


def lagrange_inversion_theorem(eq, no_of_terms = 3):
    """
    The Lagrange inversion theorem (or Lagrange inversion formula, which we abbreviate
    as LIT), also known as the Lagrange--Bürmann formula, gives the Taylor series
    expansion of the inverse function of an analytic function.

    The theorem was proved by Lagrange (1736--1813) and generalized by the German
    mathematician and teacher Hans Heinrich Bürmann ( --1817), both in the late 18th
    century.

    The Lagrange inversion formula is one of the fundamental formulas of combinatorics.
    In its simplest form it gives a formula for the power series coefficients of the
    solution f(x) of the function equation f(x) = xG(f(x)) in terms of coefficients of
    powers of G.

    Theorem: Suppose z is defined as a function of w by an equation of the form

    f(w) = z

    where f is analytic at a point a and f'(a) is not equal to 0.
    Then it is possible to invert or solve the equation for w in the form of a series:

    w = a + \sum_{n = 1}^{oo} (g_n*(z - f(a))**n)/n!

    where
                      n-1
                    d            w - a       n
    g_n =  lim   ( ------- (  -----------  )   )
          w-->a       n-1     f(w) - f(a)
                   dw

    The theorem further states that this series has a non-zero radius of convergence,
    i.e., w represents an analytic function of z in a neighbourhood of z = f(a).
    This is also called reversion of series.

    Parameters
    ==========

    eq : expression, the expansion of whose inverse is to be found.

    no_of_terms : no of terms required in the expansion; default value is 3 terms.

    Examples
    ========

    >>> from sympy.series.series import lagrange_inversion_theorem
    >>> from sympy import exp
    >>> from sympy.abc import x

    Let the equation whose inverse's expansion is to be found be y = f(x)
    then the first argument to be supplied must be f(x).
    Like here, the equation whose inverse's expansion to be found is y = (x**2)*exp(x)

    >>> eq = (x**2)*exp(x)

    If no second argument is provided, the default number of terms i.e. 3, is present
    in the expansion of the inverse of the provided equation or expression (as the first argument).

    >>> lagrange_inversion_theorem(eq)
    2*(x**2*exp(x) - E)**3*exp(-3)/27 - 7*(x**2*exp(x) - E)**2*exp(-2)/54 + (x**2*exp(x) - E)*exp(-1)/3

    Here the number of terms required in the expansion is provided.

    >>> lagrange_inversion_theorem(eq, 1)
    (x**2*exp(x) - E)*exp(-1)/3
    >>> lagrange_inversion_theorem(eq, 2)
    -7*(x**2*exp(x) - E)**2*exp(-2)/54 + (x**2*exp(x) - E)*exp(-1)/3

    >>> eq = exp(x)
    >>> lagrange_inversion_theorem(eq, 2)
    -(exp(x) - 1)**2/2 + exp(x) - 1
    >>> lagrange_inversion_theorem(eq)
    (exp(x) - 1)**3/3 - (exp(x) - 1)**2/2 + exp(x) - 1

    Note
    ====

    There is a straightforward derivation using complex analysis and contour integration;
    the complex formal power series version is a consequence of knowing the formula for polynomials,
    so the theory of analytic functions may be applied. Actually, the machinery from analytic
    function theory enters only in a formal way in this proof, in that what is really needed is
    some property of the formal residue, and a more direct formal proof is available.

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Lagrange_inversion_theorem
    .. [2] http://mathworld.wolfram.com/LagrangeInversionTheorem.html
    .. [3] http://www.cfm.brown.edu/people/dobrush/am33/Mathematica/lit.html
    """
    from sympy.abc import w, x

    try:
        eq.free_symbols
    except:
        raise AssertionError("The function provided is not analytic i.e. f' = 0. Lagrange's inversion theorem "
                             "requires that the function be analytic")
    if len(eq.free_symbols) == 2:
        if not (x in eq.free_symbols and w in eq.free_symbols):
            raise AssertionError("The function provided must be in terms of 'x' or 'w'")
    if len(eq.free_symbols) == 1:
        if not (x in eq.free_symbols or w in eq.free_symbols):
            raise AssertionError("The function provided must be in terms of 'x' or 'w'")

    n = no_of_terms
    eq = eq.subs(x, w)
    z = eq

    z_diff = diff(eq, w)
    t = 0

    while (True):
        if z_diff.subs(w, t) == 0 or not (z_diff.subs(w, t)).is_real:
            t = t + 1
            continue
        else:
            break

    a = t
    f_a = z.subs(w, a)
    sum = 0

    while (n > 0):
        l = limit(diff(((w - a) / (z - f_a)) ** n, w, n - 1), w, a)
        sum = sum + (((z - f_a) ** n) / factorial(n)) * l
        n = n - 1

    sum = sum.subs(w, x)
    return sum
