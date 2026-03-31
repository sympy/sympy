from __future__ import annotations
from sympy.core.sympify import sympify

from sympy.series.limits import limit
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.function import diff

from sympy import Dummy, Add

def series(expr, x=None, x0=0, n=6, dir="+"):
    """Series expansion of expr around point `x = x0`.

    Parameters
    ==========

    expr : Expression
           The expression whose series is to be expanded.

    x : Symbol
        It is the variable of the expression to be calculated.

    x0 : Value
         The value around which ``x`` is calculated. Can be any value
         from ``-oo`` to ``oo``.

    n : Value
        The number of terms upto which the series is to be expanded.

    dir : String, optional
          The series-expansion can be bi-directional. If ``dir="+"``,
          then (x->x0+). If ``dir="-"``, then (x->x0-). For infinite
          ``x0`` (``oo`` or ``-oo``), the ``dir`` argument is determined
          from the direction of the infinity (i.e., ``dir="-"`` for
          ``oo``).

    Examples
    ========

    >>> from sympy import series, tan, oo
    >>> from sympy.abc import x
    >>> f = tan(x)
    >>> series(f, x, 2, 6, "+")
    tan(2) + (1 + tan(2)**2)*(x - 2) + (x - 2)**2*(tan(2)**3 + tan(2)) +
    (x - 2)**3*(1/3 + 4*tan(2)**2/3 + tan(2)**4) + (x - 2)**4*(tan(2)**5 +
    5*tan(2)**3/3 + 2*tan(2)/3) + (x - 2)**5*(2/15 + 17*tan(2)**2/15 +
    2*tan(2)**4 + tan(2)**6) + O((x - 2)**6, (x, 2))

    >>> series(f, x, 2, 3, "-")
    tan(2) + (2 - x)*(-tan(2)**2 - 1) + (2 - x)**2*(tan(2)**3 + tan(2))
    + O((x - 2)**3, (x, 2))

    >>> series(f, x, 2, oo, "+")
    Traceback (most recent call last):
    ...
    TypeError: 'Infinity' object cannot be interpreted as an integer

    Returns
    =======

    Expr
        Series expansion of the expression about x0

    See Also
    ========

    sympy.core.expr.Expr.series: See the docstring of Expr.series() for complete details of this wrapper.
    """
    expr = sympify(expr)
    return expr.series(x, x0, n, dir)


def lagrange_inversion(expr, x, x0=0, n=3, dir="+"):
    r"""
    The Lagrange inversion theorem (or Lagrange inversion formula, which we abbreviate
    as LIT), also known as the Lagrange--Burmann formula, gives the Taylor series
    expansion of the inverse function of an analytic function.

    The theorem was proved by Lagrange (1736--1813) and generalized by the German
    mathematician and teacher Hans Heinrich Burmann ( --1817), both in the late 18th
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

    expr : expression, the expansion of whose inverse is to be found.

    x : variable, in terms of which equation is given, and in terms of which, the
    expansion will be provided.

    x0 : point on x-axis in the neighbourhood of which, the inverse expansion is to be found.

    n : no of terms required in the expansion; default value is 3 terms.

    dir : direction to approach the point from.

    Examples
    ========

    >>> from sympy.series.series import lagrange_inversion
    >>> from sympy import exp
    >>> from sympy.abc import x

    Let the equation whose inverse's expansion is to be found be y = f(x)
    then the first argument to be supplied must be f(x).
    Like here, the equation whose inverse's expansion to be found is y = (x**2)*exp(x)

    >>> eq = (x**2)*exp(x)

    Here the number of terms required in the expansion is provided.

    >>> lagrange_inversion(eq, x, 1, 4)
    2*(x**2*exp(x) - E)**3*exp(-3)/27 - 7*(x**2*exp(x) - E)**2*exp(-2)/54 + (x**2*exp(x) - E)*exp(-1)/3 + 1

    Here the number of terms required is, by default, 3.

    >>> lagrange_inversion(eq, x, 1)
    -7*(x**2*exp(x) - E)**2*exp(-2)/54 + (x**2*exp(x) - E)*exp(-1)/3 + 1

    Here the point with respect to which the expansion is found is changed from x=1 to x=2.

    >>> lagrange_inversion(eq, x, 2)
    -7*(x**2*exp(x) - 4*exp(2))**2*exp(-4)/512 + (x**2*exp(x) - 4*exp(2))*exp(-2)/8 + 2

    >>> eq = exp(x)
    >>> lagrange_inversion(eq, x, 2)
    -(exp(x) - exp(2))**2*exp(-4)/2 + (exp(x) - exp(2))*exp(-2) + 2

    Here the number of terms and the point with respect to which the inverse expansion will be calculated are both
    not provided explicitly.

    >>> lagrange_inversion(eq, x)
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

    .. [1] https://en.wikipedia.org/wiki/lagrange_inversion
    .. [2] http://mathworld.wolfram.com/LagrangeInversionTheorem.html
    .. [3] http://www.cfm.brown.edu/people/dobrush/am33/Mathematica/lit.html
    """

    w = sympify(expr)

    if not w.free_symbols:
        raise ValueError("The expression provided has no free symbols.")

    x = sympify(x)

    if x not in w.free_symbols:
        raise ValueError("The variable provided does not match the expression variable.")

    x0 = sympify(x0)

    if diff(w, x).subs(x, x0) == 0:
        raise ValueError("f'(a) = 0. Lagrange's inversion theorem requires "
                         "that the derivative be non-zero at point 'x0'.")

    if n <= 0:
        return x0

    # Precomputing
    fx0 = w.subs(x, x0)
    t = Dummy('t')
    w_t = w.subs(x, t)

    # Construct the series using a single generator expression
    # Each term k corresponds to (1/k!) * [ d^(k-1)/dt^(k-1) ((t-x0)/(f(t)-f(x0)))^k ]_{t=x0} * (f(x)-f(x0))^k
    return Add(x0, *(
        (limit(diff(((t - x0) / (w_t - fx0))**k, t, k - 1), t, x0, dir=dir) / factorial(k)) * (w - fx0)**k
        for k in range(1, n)
    ), evaluate=True)
