from sympy.core.sympify import sympify
from sympy import symbols


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
          then (x->x0+). If ``dir="-", then (x->x0-). For infinite
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


def inversion_series(f, x, a=0, n=3):
    r"""
    The Lagrange inversion theorem (or Lagrange inversion formula, which we abbreviate
    as LIT), also known as the Lagrange--Bürmann formula, gives the Taylor series
    expansion of the inverse function of an analytic function.
    It finds important applications in combinatorics and calculus.

    The theorem was proved by Lagrange (1736--1813) and generalized by the German
    mathematician and teacher Hans Heinrich Bürmann ( --1817), both in the late 18th
    century.

    Theorem: Suppose x is defined as a function of y by an equation of the form

    f(x) = y

    where f is analytic at a point a and f'(a) is not equal to 0.
    Then it is possible to invert or solve the equation for x in the form of a series:

    .. math::
        g(y) = a + \sum_{n = 1}^{oo} \frac{g_n (y - f(a))^n}{n!}

    where

    .. math::
        g_n = \lim_{y \to a} \left[ \frac{d^{n-1}}{dx^{n-1}}
        \left( \frac{x-a}{f(x)-f(a)} \right) ^n \right]

    The theorem further states that this series has a non-zero radius of convergence,
    i.e., x represents an analytic function of y in a neighbourhood of y = f(a).

    Parameters
    ==========

    f : expression, the expansion of whose inverse is to be found.

    x : variable, in terms of which equation is given, and in terms of which, the
    expansion will be provided.

    a : point on x-axis in the neighbourhood of which, the inverse expansion is to be found.

    n : no of terms required in the expansion; default value is 3 terms.

    Examples
    ========

    >>> from sympy.series.series import inversion_series
    >>> from sympy import exp
    >>> from sympy.abc import x

    Let the equation whose inverse's expansion is to be found be y = f(x)
    then the first argument to be supplied must be f(x).

    >>> f = exp(x)

    In this case, the inverse series is that of the logarithmic function:

    >>> inversion_series(f,x,0,4)
    x + (x-1)**3/3 - (x-1)**2/2 - 1

    Functions with vanishing derivative cannot be inverted with this method
    >>> f = x**2
    >>> inversion_series(f,x,0,2)
    ValueError: The inversion formula requires non-zero derivative.

    Note
    ====

    There is a straightforward derivation using complex analysis and contour integration;
    the complex formal power series version is a consequence of knowing the formula for polynomials,
    so the theory of analytic functions may be applied. Actually, the machinery from analytic
    function theory enters only in a formal way in this proof, in that what is really needed is
    some property of the formal residue, and a more direct formal proof is available.

    The implementation uses a quasi-newton method, iterating x:=(f(x)-y)/f'(a)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Lagrange_inversion_theorem
    .. [2] http://mathworld.wolfram.com/LagrangeInversionTheorem.html
    .. [3] http://www.cfm.brown.edu/people/dobrush/am33/Mathematica/lit.html
       [4] Abramowitz, M. and Stegun, I.A., 1964. "Handbook of Mathematical Functions", Eq. 3.6.6
    """

    if f.diff(x,1).subs({x:a}) == 0:
        raise ValueError("The inversion formula requires non-zero derivative.")

    # solve using quasi-newton method
    fa = f.subs({x:a})
    f = f.series(x,a,n).removeO()
    dfa = f.diff(x,1).subs({x:a})

    y = symbols('y')
    x_series = a + (y-fa)/dfa

    for i in range(2,n):
        step = f.subs({x:x_series}).series(y,fa,i+1).coeff(y-fa,i)
        x_series = x_series-step*(y-fa)**i/dfa

    x_series = x_series.subs({y:x})

    return x_series
