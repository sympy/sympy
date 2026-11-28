from __future__ import annotations
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.symbol import Dummy
from sympy.functions import SingularityFunction, DiracDelta
from sympy.integrals import integrate
from sympy.polys.polytools import Poly


def singularityintegrate(f, x):
    """
    This function handles the indefinite integrations of Singularity functions.
    The ``integrate`` function calls this function internally whenever an
    instance of SingularityFunction is passed as argument.

    Explanation
    ===========

    The idea for integration is the following:

    - If we are dealing with a SingularityFunction expression,
      i.e. ``SingularityFunction(x, a, n)``, we just return
      ``SingularityFunction(x, a, n + 1)/(n + 1)`` if ``n >= 0`` and
      ``SingularityFunction(x, a, n + 1)`` if ``n < 0``.

    - If the node is a multiplication or power node having a
      SingularityFunction term we rewrite the whole expression in terms of
      Heaviside and DiracDelta and then integrate the output. Lastly, we
      rewrite the output of integration back in terms of SingularityFunction.

    - If none of the above case arises, we return None.

    Examples
    ========

    >>> from sympy.integrals.singularityfunctions import singularityintegrate
    >>> from sympy import SingularityFunction, symbols, Function
    >>> x, a, n, y = symbols('x a n y')
    >>> f = Function('f')
    >>> singularityintegrate(SingularityFunction(x, a, 3), x)
    SingularityFunction(x, a, 4)/4
    >>> singularityintegrate(5*SingularityFunction(x, 5, -2), x)
    5*SingularityFunction(x, 5, -1)
    >>> singularityintegrate(6*SingularityFunction(x, 5, -1), x)
    6*SingularityFunction(x, 5, 0)
    >>> singularityintegrate(x*SingularityFunction(x, 0, -1), x)
    0
    >>> singularityintegrate(SingularityFunction(x, 1, -1) * f(x), x)
    f(1)*SingularityFunction(x, 1, 0)

    A polynomial times a SingularityFunction with a nonnegative exponent is
    integrated in terms of SingularityFunctions:

    >>> singularityintegrate((x - 4)*SingularityFunction(x, 4, 0), x)
    SingularityFunction(x, 4, 2)/2

    """

    if not f.has(SingularityFunction):
        return None

    if isinstance(f, SingularityFunction):
        x, a, n = f.args
        if n.is_positive or n.is_zero:
            return SingularityFunction(x, a, n + 1)/(n + 1)
        elif n in (-1, -2, -3, -4):
            return SingularityFunction(x, a, n + 1)

    if f.is_Mul:
        rv = _poly_singularityintegrate(f, x)
        if rv is not None:
            return rv

    if f.is_Mul or f.is_Pow:

        expr = f.rewrite(DiracDelta)
        expr = integrate(expr, x)
        return expr.rewrite(SingularityFunction)
    return None


def _poly_singularityintegrate(f, x):
    """
    Integrate ``p(x)*SingularityFunction(x, a, n)`` where ``p`` is a
    polynomial in ``x`` and ``n`` is nonnegative.

    Writing ``p(x)`` as a polynomial in ``x - a``, i.e.
    ``p(x) = Sum(c_k*(x - a)**k)``, gives
    ``p(x)*<x - a>**n = Sum(c_k*<x - a>**(n + k))`` because
    ``<x - a>**n = (x - a)**n*Heaviside(x - a)`` for ``n >= 0``. Each term can
    then be integrated directly. Return None if ``f`` is not of this form.
    """
    sfs = [arg for arg in f.args if isinstance(arg, SingularityFunction)]
    if len(sfs) != 1:
        return None
    sf = sfs[0]
    var, a, n = sf.args
    if var != x or a.has(x) or not n.is_nonnegative:
        return None
    rest = Mul(*[arg for arg in f.args if arg is not sf])
    if rest.has(SingularityFunction) or not rest.is_polynomial(x):
        return None
    t = Dummy('t')
    p = Poly(rest.subs(x, t + a), t)
    return Add(*[c*SingularityFunction(x, a, n + k + 1)/(n + k + 1)
                 for (k,), c in p.terms()])
