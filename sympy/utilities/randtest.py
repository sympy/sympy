""" Helpers for randomised testing """

from sympy import I
from random import uniform

def random_complex_number():
    """
    Return a random complex number.

    To reduce chance of hitting branch cuts or anything, we guarantee
    -1 <= Im z <= 1, 2 <= Re z <= 3
    """
    return uniform(2, 3) + I*uniform(-1, 1)

def test_numerically(f, g, z, tol=1.0e-6):
    """
    Test numerically that f and g agree when evaluated in the argument z.

    Examples:

    >>> from sympy import sin, cos, S
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import test_numerically as tn
    >>> tn(sin(x)**2 + cos(x)**2, S(1), x)
    True
    """
    z0 = random_complex_number()
    return abs((f.subs(z, z0) - g.subs(z, z0)).n()) <= tol

def test_derivative_numerically(f, z, tol=1.0e-6):
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
    z0 = random_complex_number()
    f1 = f.diff(z).subs(z, z0)
    f2 = Derivative(f, z).doit_numerically(z0)
    return abs((f1 - f2).n()) <= tol
