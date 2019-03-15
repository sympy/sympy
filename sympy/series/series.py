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


def lagrange_inversion_theorem(eq, x, a=0, n=3):
    r"""
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

    .. math::
        g(z) = a + \sum_{n = 1}^{oo} \frac{g_n (z - f(a))^n}{n!}

    where

    .. math::
        g_n = \lim_{w \to a} \left[ \frac{d^{n-1}}{dw^{n-1}}
        \left( \frac{w-a}{f(w)-f(a)} \right) ^n \right]

    The theorem further states that this series has a non-zero radius of convergence,
    i.e., w represents an analytic function of z in a neighbourhood of z = f(a).
    This is also called reversion of series.

    Parameters
    ==========

    eq : expression, the expansion of whose inverse is to be found.

    x : variable, in terms of which equation is given, and in terms of which, the
    expansion will be provided.

    a : point on x-axis in the neighbourhood of which, the inverse expansion is to be found.

    n : no of terms required in the expansion; default value is 3 terms.

    Examples
    ========

    >>> from sympy.series.series import lagrange_inversion_theorem
    >>> from sympy import exp
    >>> from sympy.abc import x

    Let the equation whose inverse's expansion is to be found be y = f(x)
    then the first argument to be supplied must be f(x).
    Like here, the equation whose inverse's expansion to be found is y = (x**2)*exp(x)

    >>> eq = (x**2)*exp(x)

    Here the number of terms required in the expansion is provided.

    >>> lagrange_inversion_theorem(eq, x, 1, 4)
    2*(x**2*exp(x) - E)**3*exp(-3)/27 - 7*(x**2*exp(x) - E)**2*exp(-2)/54 + (x**2*exp(x) - E)*exp(-1)/3 + 1

    Here the number of terms required is, by default, 3.

    >>> lagrange_inversion_theorem(eq, x, 1)
    -7*(x**2*exp(x) - E)**2*exp(-2)/54 + (x**2*exp(x) - E)*exp(-1)/3 + 1

    Here the point with respect to which the expansion is found is changed from x=1 to x=2.

    >>> lagrange_inversion_theorem(eq, x, 2)
    -7*(x**2*exp(x) - 4*exp(2))**2*exp(-4)/512 + (x**2*exp(x) - 4*exp(2))*exp(-2)/8 + 2

    >>> eq = exp(x)
    >>> lagrange_inversion_theorem(eq, x, 2)
    -(exp(x) - exp(2))**2*exp(-4)/2 + (exp(x) - exp(2))*exp(-2) + 2

    Here the number of terms and the point with respect to which the inverse expansion will be calculated are both
    not provided explicitly.

    >>> lagrange_inversion_theorem(eq, x)
    -(exp(x) - 1)**2/2 + exp(x) - 1

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

    try:
        not eq.free_symbols
    except AttributeError:
        raise ValueError("The function provided is not analytic i.e. f' = 0. Lagrange's inversion theorem requires "
                         "that the function be analytic")

    if len(eq.free_symbols) == 1:
        if x not in eq.free_symbols:
            raise ValueError("The function variable provided does not match the equation variable")

    w = x
    eq = eq.subs(x, w)
    z = eq

    z_diff = diff(eq, w)

    if z_diff.subs(w, a) == 0 or not (z_diff.subs(w, a)).is_real:
        raise ValueError("The first derivative of the function is zero at the point 'a' provided i.e. f'(a) = 0. "
                         "Lagrange's inversion theorem requires that the derivative of the function be non zero at "
                         "the point 'a'.")

    f_a = z.subs(w, a)
    result = a

    if n == 1:
        return result

    a_1 = ((w - a)/(z - f_a))
    l = limit(a_1, w, a)
    result = result + l*(z - f_a)

    if n == 2:
        return result

    for i in range(2, n):
        l = limit(diff(((w - a) / (z - f_a)) ** i, w, i - 1), w, a)
        result = result + (((z - f_a) ** i) / factorial(i)) * l

    result = result.subs(w, x)
    return result
