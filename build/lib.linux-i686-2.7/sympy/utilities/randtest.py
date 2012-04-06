""" Helpers for randomized testing """

from sympy import I, nsimplify, S, Tuple, Dummy
from random import uniform

def random_complex_number(a=2, b=-1, c=3, d=1, rational=False):
    """
    Return a random complex number.

    To reduce chance of hitting branch cuts or anything, we guarantee
    b <= Im z <= d, a <= Re z <= c
    """
    A, B = uniform(a, c), uniform(b, d)
    if not rational:
        return A + I*B
    return nsimplify(A, rational=True) + I*nsimplify(B, rational=True)


def comp(z1, z2, tol):
    """Return a bool indicating whether the error between z1 and z2 is <= tol.

    If z2 is non-zero and ``|z1| > 1`` the error is normalized by ``|z1|``, so
    if you want the absolute error, call this as ``comp(z1 - z2, 0, tol)``.
    """
    if not z1:
        z1, z2 = z2, z1
    if not z1:
        return True
    diff = abs(z1 - z2)
    az1 = abs(z1)
    if z2 and az1 > 1:
        return diff/az1 <= tol
    else:
        return diff <= tol

def test_numerically(f, g, z=None, tol=1.0e-6, a=2, b=-1, c=3, d=1):
    """
    Test numerically that f and g agree when evaluated in the argument z.

    If z is None, all symbols will be tested. This routine does not test
    whether there are Floats present with precision higher than 15 digits
    so if there are, your results may not be what you expect due to round-
    off errors.

    Examples
    ========

    >>> from sympy import sin, cos, S
    >>> from sympy.abc import x
    >>> from sympy.utilities.randtest import test_numerically as tn
    >>> tn(sin(x)**2 + cos(x)**2, 1, x)
    True
    """
    f, g, z = Tuple(f, g, z)
    z = [z] if z else (f.free_symbols | g.free_symbols)
    reps = zip(z, [random_complex_number(a, b, c, d) for zi in z])
    z1 = f.subs(reps).n()
    z2 = g.subs(reps).n()
    return comp(z1, z2, tol)

def test_derivative_numerically(f, z, tol=1.0e-6, a=2, b=-1, c=3, d=1):
    """
    Test numerically that the symbolically computed derivative of f
    with respect to z is correct.

    This routine does not test whether there are Floats present with
    precision higher than 15 digits so if there are, your results may
    not be what you expect due to round-off errors.

    Examples
    ========

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
