""" Helpers for randomised testing """

from sympy import I
from random import uniform

def random_complex_number(a=2, b=-1, c=3, d=1):
    """
    Return a random complex number.

    To reduce chance of hitting branch cuts or anything, we guarantee
    b <= Im z <= d, a <= Re z <= c
    """
    return uniform(a, c) + I*uniform(b, d)

def comp(z1, z2, tol):
    diff = abs(z1 - z2)
    if abs(z1) > 1:
        return diff/abs(z1) <= tol
    else:
        return diff <= tol

def test_numerically(f, g, z, tol=1.0e-6, a=2, b=-1, c=3, d=1):
    """
    Test numerically that f and g agree when evaluated in the argument z.

    Examples:

    >>> from sympy import sin, cos, S
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import test_numerically as tn
    >>> tn(sin(x)**2 + cos(x)**2, S(1), x)
    True
    """
    z0 = random_complex_number(a, b, c, d)
    z1 = f.subs(z, z0).n()
    z2 = g.subs(z, z0).n()
    return comp(z1, z2, tol)

def test_derivative_numerically(f, z, tol=1.0e-6, a=2, b=-1, c=3, d=1):
    """
    Test numerically that the symbolically computed derivative of f
    with respect to z is correct.

    Examples:
    >>> from sympy import sin, cos
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import test_derivative_numerically as td
    >>> td(sin(x), x)
    True
    """
    from sympy.core.function import Derivative
    z0 = random_complex_number(a, b, c, d)
    f1 = f.diff(z).subs(z, z0)
    f2 = Derivative(f, z).doit_numerically(z0)
    return comp(f1.n(), f2.n(), tol)
