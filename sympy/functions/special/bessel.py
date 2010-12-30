from math import pi

from sympy.functions.elementary.trigonometric import sin, cos
from sympy.parsing.sympify import sympify

def fn(n, z):
    """
    Coefficients for the spherical Bessel functions.

    Those are only needed in the jn() function.

    The coefficients are calculated from:

    fn(0, z) = 1/z
    fn(1, z) = 1/z**2
    fn(n-1, z) + fn(n+1, z) == (2*n+1)/z * fn(n, z)

    Examples:

    >>> from sympy.functions.special.bessel import fn
    >>> from sympy import Symbol
    >>> z = Symbol("z")
    >>> fn(1, z)
    z**(-2)
    >>> fn(2, z)
    -1/z + 3/z**3
    >>> fn(3, z)
    15/z**4 - 6/z**2
    >>> fn(4, z)
    1/z + 105/z**5 - 45/z**3

    """
    n = sympify(n)
    if not n.is_Integer:
        raise TypeError("'n' must be an Integer")
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

        >>> from sympy import Symbol, jn, sin, cos
        >>> z = Symbol("z")
        >>> print jn(0, z)
        sin(z)/z
        >>> jn(1, z) == sin(z)/z**2 - cos(z)/z
        True
        >>> jn(3, z) ==(1/z - 15/z**3)*cos(z) + (15/z**4 - 6/z**2)*sin(z)
        True

    The spherical Bessel functions are calculated using the formula:

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

        >>> from sympy import Symbol, yn, sin, cos
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

def jn_zeros(n, k, method="sympy"):
    """
    Zeros of the spherical Bessel function of the first kind.

    This returns an array of zeros of jn up to the k-th zero.

    method = "sympy": uses the SymPy's jn and findroot to find all roots
    method = "scipy": uses the SciPy's sph_jn and newton to find all roots,
            which if faster than method="sympy", but it requires SciPy and only
            works with low precision floating point numbers

    Examples:

        >>> from sympy.mpmath import nprint
        >>> from sympy import jn_zeros
        >>> nprint(jn_zeros(2, 4))
        [5.76345919689, 9.09501133048, 12.3229409706, 15.5146030109]

    """
    if method == "sympy":
        from sympy.mpmath import findroot
        f = lambda x: jn(n, x).n()
    elif method == "scipy":
        from scipy.special import sph_jn
        from scipy.optimize import newton
        f  = lambda x: sph_jn(n, x)[0][-1]
    elif method == 'mpmath':
        # this needs a recent version of mpmath, newer than in sympy
        from mpmath import besseljzero
        return [besseljzero(n + 0.5, k) for k in xrange(1, k + 1)]
    else:
        raise NotImplementedError("Unknown method.")
    def solver(f, x):
        if method == "sympy":
            # findroot(solver="newton") or findroot(solver="secant") can't find
            # the root within the given tolerance. So we use solver="muller",
            # which converges towards complex roots (even for real starting
            # points), and so we need to chop all complex parts (that are small
            # anyway). Also we need to set the tolerance, as it sometimes fail
            # without it.
            def f_real(x):
                return f(complex(x).real)
            root = findroot(f_real, x, solver="muller", tol=1e-9)
            root = complex(root).real
        elif method == "scipy":
            root = newton(f, x)
        else:
            raise NotImplementedError("Unknown method.")
        return root

    # we need to approximate the position of the first root:
    root = n+pi
    # determine the first root exactly:
    root = solver(f, root)
    roots = [root]
    for i in range(k-1):
        # estimate the position of the next root using the last root + pi:
        root = solver(f, root+pi)
        roots.append(root)
    return roots
