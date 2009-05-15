from sympy.core import sympify
from sympy.functions.elementary.trigonometric import sin, cos

def fn(n, z):
    """
    Coefficients for the spherical Bessel functions.

    Those are only needed in the jn() function.

    The coefficients are calculated from:

    fn(0, z) = 1/z
    fn(1, z) = 1/z**2
    fn(n-1, z) + fn(n+1, z) == (2*n+1)/z * fn(n, z)
    """
    if n == 0:
        return 1/z
    elif n == 1:
        return 1/z**2
    elif n > 1:
        return ((2*n-1)/z * fn(n-1, z)).expand() - fn(n-2, z)
    elif n < 0:
        return ((2*n+3)/z * fn(n+1, z)).expand() - fn(n+2, z)


def jn(n, z):
    """
    Spherical Bessel function of the first kind.

    Examples:

        >>> from sympy import Symbol
        >>> z = Symbol("z")
        >>> print jn(0, z)
        sin(z)/z
        >>> jn(1, z) == sin(z)/z**2 - cos(z)/z
        True
        >>> jn(3, z) ==(1/z - 15/z**3)*cos(z) + (15/z**4 - 6/z**2)*sin(z)
        True

    The spherical Bessel functions are calculated using the formula[0]:

    jn(n, z) == fn(n, z) * sin(z) + (-1)**(n+1) * fn(-n-1, z) * cos(z)

    where fn(n, z) are the coefficients, see fn()'s sourcecode for more
    information.
    """

    n = sympify(n)
    z = sympify(z)
    return fn(n, z) * sin(z) + (-1)**(n+1) * fn(-n-1, z) * cos(z)

def yn(n, z):
    """
    Spherical Bessel function of the second kind.

    Examples:

        >>> from sympy import Symbol
        >>> z = Symbol("z")
        >>> print yn(0, z)
        -cos(z)/z
        >>> yn(1, z) == -cos(z)/z**2-sin(z)/z
        True

    yn is calculated using the formula:

    yn(n, z) == (-1)**(n+1) * jn(-n-1, z)
    """

    n = sympify(n)
    z = sympify(z)
    return (-1)**(n+1) * jn(-n-1, z)
