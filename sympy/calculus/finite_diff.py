"""
Finite difference weights
=========================

This module implements an algorithm for efficient generation of finite
difference weights for ordinary differentials of functions for
derivatives from 0 (interpolation) up to arbitrary order.

The core algorithm is provided in the finite difference weight generating
method (finite_diff_weights), and a convenience function for estimating
a derivative (or interpolate) directly from a series of points
is also provided (apply_finite_diff).

"""

from sympy import S


def finite_diff_weights(order, x_list, x0):
    """
    Calculates the finite difference weights for an arbitrarily
    spaced one-dimensional grid (x_list) for derivatives at 'x0'
    of order 0, 1, ..., up to 'order' using a recursive formula.

    Parameters
    ==========

    order : int
        Up to what derivative order weights should be calculated.
        0 corresponds to interpolation.
    x_list: sequence
        Strictly monotonically increasing sequence of values for
        the independent variable.
    x0: Number or Symbol
        At what value of the independent variable the finite difference
        weights should be generated.

    Returns
    =======

    list
        A list of sublists, each corresponding to coefficients for
        increasing derivative order, and each containing lists of
        coefficients for increasing accuracy.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.calculus import finite_diff_weights
    >>> finite_diff_weights(1, [-S(1)/2, S(1)/2, S(3)/2, S(5)/2], 0)
    [[[1, 0, 0, 0],
      [1/2, 1/2, 0, 0],
      [3/8, 3/4, -1/8, 0],
      [5/16, 15/16, -5/16, 1/16]],
     [[0, 0, 0, 0], [-1, 1, 0, 0], [-1, 1, 0, 0], [-23/24, 7/8, 1/8, -1/24]]]

    the result is two subslists, the first is for the 0:th derivative
    (interpolation) and the second for the first derivative (we gave
    1 as the parameter of order so this is why we get no list for
    a higher order derivative). Each sublist contains the most accurate
    formula in the end (all points used).

    Beware of the offset in the lower accuracy formulae when looking at a
    centered difference:

    >>> from sympy import S
    >>> from sympy.calculus import finite_diff_weights
    >>> finite_diff_weights(1, [-S(5)/2, -S(3)/2, -S(1)/2, S(1)/2,\
    S(3)/2, S(5)/2], 0)
    [[[1, 0, 0, 0, 0, 0],
      [-3/2, 5/2, 0, 0, 0, 0],
      [3/8, -5/4, 15/8, 0, 0, 0],
      [1/16, -5/16, 15/16, 5/16, 0, 0],
      [3/128, -5/32, 45/64, 15/32, -5/128, 0],
      [3/256, -25/256, 75/128, 75/128, -25/256, 3/256]],
     [[0, 0, 0, 0, 0, 0],
      [-1, 1, 0, 0, 0, 0],
      [1, -3, 2, 0, 0, 0],
      [1/24, -1/8, -7/8, 23/24, 0, 0],
      [0, 1/24, -9/8, 9/8, -1/24, 0],
      [-3/640, 25/384, -75/64, 75/64, -25/384, 3/640]]]


    The capability to generate weights at arbitrary points can be
    used e.g. to minimize Runge's phenomenon by using Chebyshev nodes:

    >>> from sympy import cos, symbols, pi, simplify
    >>> from sympy.calculus import finite_diff_weights
    >>> N, (h, x) = 4, symbols('h x')
    >>> x_list = [x+h*cos(i*pi/(N)) for i in range(N,-1,-1)] # chebyshev nodes
    >>> print(x_list)
    [-h + x, -sqrt(2)*h/2 + x, x, sqrt(2)*h/2 + x, h + x]
    >>> mycoeffs = finite_diff_weights(1, x_list, 0)[1][4]
    >>> [simplify(c) for c in  mycoeffs]
    [(h**3/2 + h**2*x - 3*h*x**2 - 4*x**3)/h**4,\
    (-sqrt(2)*h**3 - 4*h**2*x + 3*sqrt(2)*h*x**2 + 8*x**3)/h**4,\
    6*x/h**2 - 8*x**3/h**4,
    (sqrt(2)*h**3 - 4*h**2*x - 3*sqrt(2)*h*x**2 + 8*x**3)/h**4,\
    (-h**3/2 + h**2*x + 3*h*x**2 - 4*x**3)/h**4]

    Notes
    =====

    If weights for a finite difference approximation
    of the 3rd order derivative is wanted, weights for 0th, 1st
    and 2nd order are calculated "for free", so are formulae using
    fewer and fewer of the paramters. This is something one can
    take advantage of to save computational cost.

    See also
    ========

    sympy.calculus.finite_diff.apply_finite_diff


    References
    ==========

    .. [1] Generation of Finite Difference Formulas on Arbitrarily Spaced
            Grids, Bengt Fornberg; Mathematics of compuation; 51; 184;
            (1988); 699-706; doi:10.1090/S0025-5718-1988-0935077-0

    """
    # The notation below colosely corresponds that used in the paper.
    if order < 0:
        raise ValueError("Negtive derivative order illegal.")
    if int(order) != order:
        raise ValueError("Non-integer order illegal")
    M = order
    N = len(x_list) - 1
    delta = [[[0 for nu in range(N+1)] for n in range(N+1)] for
             m in range(M+1)]
    delta[0][0][0] = S(1)
    c1 = S(1)
    for n in range(1, N+1):
        c2 = S(1)
        for nu in range(0, n):
            c3 = x_list[n]-x_list[nu]
            c2 = c2 * c3
            if n <= M:
                delta[n][n-1][nu] = 0
            for m in range(0, min(n, M)+1):
                delta[m][n][nu] = (x_list[n]-x0)*delta[m][n-1][nu] -\
                    m*delta[m-1][n-1][nu]
                delta[m][n][nu] /= c3
        for m in range(0, min(n, M)+1):
            delta[m][n][n] = c1/c2*(m*delta[m-1][n-1][n-1] -
                                    (x_list[n-1]-x0)*delta[m][n-1][n-1])
        c1 = c2
    return delta


def apply_finite_diff(order, x_list, y_list, x0):
    """
    Calculates the finite difference approxiamtion of
    the derivative of requested order at x0 from points
    provided in x_list and y_list.

    Parameters
    ==========

    order: int
        order of derivative to approximate. 0 corresponds to interpolation.
    x_list: sequence
        Strictly monotonically increasing sequence of values for
        the independent variable.
    y_list: sequence
        The function value at corresponding values for the independent
        variable in x_list.
    x0: Number or Symbol
        At what value of the independent variable the derivative should be
        evaluated.

    Returns
    =======

    sympy.core.add.Add or sympy.core.numbers.Number
        The finite difference expression approximating the requested
        derivative order at x0.

    Examples
    ========

    >>> from sympy.calculus import apply_finite_diff
    >>> cube = lambda arg: (1.0*arg)**3
    >>> xlist = range(-3,3+1)
    >>> apply_finite_diff(2, xlist, map(cube, xlist), 2) - 12 # doctest: +SKIP
    -3.55271367880050e-15

    we see that the example above only contain rounding errors.
    apply_finite_diff can also be used on more abstract objects:

    >>> from sympy import IndexedBase, Idx
    >>> from sympy.calculus import apply_finite_diff
    >>> x, y = map(IndexedBase, 'xy')
    >>> i = Idx('i')
    >>> x_list, y_list = zip(*[(x[i+j], y[i+j]) for j in range(-1,2)])
    >>> apply_finite_diff(1, x_list, y_list, x[i])
    (-1 + (x[i + 1] - x[i])/(-x[i - 1] + x[i]))*y[i]/(x[i + 1] - x[i]) + \
    (-x[i - 1] + x[i])*y[i + 1]/((-x[i - 1] + x[i + 1])*(x[i + 1] - x[i])) -\
    (x[i + 1] - x[i])*y[i - 1]/((-x[i - 1] + x[i + 1])*(-x[i - 1] + x[i]))


    Notes
    =====

    Order = 0 corresponds to interpolation.
    Only supply so many points you think makes sense
    to around x0 when extracting the derivative (the function
    need to be well behaved within that region). Also beware
    of Runge's phenomenon.

    See also
    ========

    sympy.calculus.finite_diff.finite_diff_weights

    References
    ==========

    Fortran 90 implementation with Python interface for numerics: finitediff_

    .. _finitediff: https://github.com/bjodah/finitediff

    """

    # In the original paper the following holds for the notation:
    # M = order
    # N = len(x_list) - 1

    N = len(x_list) - 1
    if len(x_list) != len(y_list):
        raise ValueError("x_list and y_list not equal in length.")

    delta = finite_diff_weights(order, x_list, x0)

    derivative = 0
    for nu in range(0, len(x_list)):
        derivative += delta[order][N][nu]*y_list[nu]
    return derivative
