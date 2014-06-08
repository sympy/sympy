#!/usr/bin/env python

"""
This example illustrates the Gibbs phenomenon.

It also calculates the Wilbraham-Gibbs constant by two approaches:

1) calculating the fourier series of the step function and determining the
first maximum.
2) evaluating the integral for si(pi).

See:
 * http://en.wikipedia.org/wiki/Gibbs_phenomena
"""

from sympy import var, sqrt, integrate, conjugate, seterr, Abs, pprint, I, pi,\
    sin, cos, sign, Plot, lambdify, Integral, S

# seterr(True)

x = var("x", real=True)


def l2_norm(f, lim):
    """
    Calculates L2 norm of the function "f", over the domain lim=(x, a, b).

    x ...... the independent variable in f over which to integrate
    a, b ... the limits of the interval

    Example:

    >>> from sympy import Symbol
    >>> from gibbs_phenomenon import l2_norm
    >>> x = Symbol('x', real=True)
    >>> l2_norm(1, (x, -1, 1))
    sqrt(2)
    >>> l2_norm(x, (x, -1, 1))
    sqrt(6)/3

    """
    return sqrt(integrate(Abs(f)**2, lim))


def l2_inner_product(a, b, lim):
    """
    Calculates the L2 inner product (a, b) over the domain lim.
    """
    return integrate(conjugate(a)*b, lim)


def l2_projection(f, basis, lim):
    """
    L2 projects the function f on the basis over the domain lim.
    """
    r = 0
    for b in basis:
        r += l2_inner_product(f, b, lim) * b
    return r


def l2_gram_schmidt(list, lim):
    """
    Orthonormalizes the "list" of functions using the Gram-Schmidt process.

    Example:

    >>> from sympy import Symbol
    >>> from gibbs_phenomenon import l2_gram_schmidt

    >>> x = Symbol('x', real=True)    # perform computations over reals to save time
    >>> l2_gram_schmidt([1, x, x**2], (x, -1, 1))
    [sqrt(2)/2, sqrt(6)*x/2, 3*sqrt(10)*(x**2 - 1/3)/4]

    """
    r = []
    for a in list:
        if r == []:
            v = a
        else:
            v = a - l2_projection(a, r, lim)
        v_norm = l2_norm(v, lim)
        if v_norm == 0:
            raise ValueError("The sequence is not linearly independent.")
        r.append(v/v_norm)
    return r


def integ(f):
    return integrate(f, (x, -pi, 0)) + integrate(-f, (x, 0, pi))


def series(L):
    """
    Normalizes the series.
    """
    r = 0
    for b in L:
        r += integ(b)*b
    return r


def msolve(f, x):
    """
    Finds the first root of f(x) to the left of 0.

    The x0 and dx below are taylored to get the correct result for our
    particular function --- the general solver often overshoots the first
    solution.
    """
    f = lambdify(x, f)
    x0 = -0.001
    dx = 0.001
    while f(x0 - dx) * f(x0) > 0:
        x0 = x0 - dx
    x_max = x0 - dx
    x_min = x0
    assert f(x_max) > 0
    assert f(x_min) < 0
    for n in range(100):
        x0 = (x_max + x_min)/2
        if f(x0) > 0:
            x_max = x0
        else:
            x_min = x0
    return x0


def main():
    # L = l2_gram_schmidt([1, cos(x), sin(x), cos(2*x), sin(2*x)], (x, -pi, pi))
    # L = l2_gram_schmidt([1, cos(x), sin(x)], (x, -pi, pi))
    # the code below is equivalen to l2_gram_schmidt(), but faster:
    L = [1/sqrt(2)]
    for i in range(1, 100):
        L.append(cos(i*x))
        L.append(sin(i*x))
    L = [f/sqrt(pi) for f in L]

    f = series(L)
    print("Fourier series of the step function")
    pprint(f)
    # Plot(f.diff(x), [x, -5, 5, 3000])
    x0 = msolve(f.diff(x), x)

    print("x-value of the maximum:", x0)
    max = f.subs(x, x0).evalf()
    print("y-value of the maximum:", max)
    g = max*pi/2
    print("Wilbraham-Gibbs constant        :", g.evalf())
    print("Wilbraham-Gibbs constant (exact):", \
        Integral(sin(x)/x, (x, 0, pi)).evalf())

if __name__ == "__main__":
    main()
