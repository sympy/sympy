"""
Finite difference weights
=========================

This module implements an algorithm for efficient generation of finite
difference weights for ordinary differentials of functions for
derivatives from 0 (interpolation) up to arbitrary order.

The core algorithm is provided in the finite difference weight generating
function (finite_diff_weights), and two convenience functions are provided
for:

- estimating a derivative (or interpolate) directly from a series of points
    is also provided (``apply_finite_diff``).
- making a finite difference approximation of a Derivative instance
    (``as_finite_diff``).

"""

from sympy import S
from sympy.core.compatibility import iterable, range


def finite_diff_weights(order, x_list, x0=S(0)):
    """
    Calculates the finite difference weights for an arbitrarily
    spaced one-dimensional grid (x_list) for derivatives at 'x0'
    of order 0, 1, ..., up to 'order' using a recursive formula.
    Order of accuracy is at least len(x_list) - order, if x_list
    is defined accurately.


    Parameters
    ==========

    order: int
        Up to what derivative order weights should be calculated.
        0 corresponds to interpolation.
    x_list: sequence
        Sequence of (unique) values for the independent variable.
        It is usefull (but not necessary) to order x_list from
        nearest to farest from x0; see examples below.
    x0: Number or Symbol
        Root or value of the independent variable for which the finite
        difference weights should be generated. Defaults to S(0).

    Returns
    =======

    list
        A list of sublists, each corresponding to coefficients for
        increasing derivative order, and each containing lists of
        coefficients for increasing subsets of x_list.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.calculus import finite_diff_weights
    >>> res = finite_diff_weights(1, [-S(1)/2, S(1)/2, S(3)/2, S(5)/2], 0)
    >>> res
    [[[1, 0, 0, 0],
      [1/2, 1/2, 0, 0],
      [3/8, 3/4, -1/8, 0],
      [5/16, 15/16, -5/16, 1/16]],
     [[0, 0, 0, 0],
      [-1, 1, 0, 0],
      [-1, 1, 0, 0],
      [-23/24, 7/8, 1/8, -1/24]]]
    >>> res[0][-1]  # FD weights for 0th derivative, using full x_list
    [5/16, 15/16, -5/16, 1/16]
    >>> res[1][-1]  # FD weights for 1st derivative
    [-23/24, 7/8, 1/8, -1/24]
    >>> res[1][-2]  # FD weights for 1st derivative, using x_list[:-1]
    [-1, 1, 0, 0]
    >>> res[1][-1][0]  # FD weight for 1st deriv. for x_list[0]
    -23/24
    >>> res[1][-1][1]  # FD weight for 1st deriv. for x_list[1], etc.
    7/8

    Each sublist contains the most accurate formula at the end.
    Note, that in the above example res[1][1] is the same as res[1][2].
    Since res[1][2] has an order of accuracy of
    len(x_list[:3]) - order = 3 - 1 = 2, the same is true for res[1][1]!

    >>> from sympy import S
    >>> from sympy.calculus import finite_diff_weights
    >>> res = finite_diff_weights(1, [S(0), S(1), -S(1), S(2), -S(2)], 0)[1]
    >>> res
    [[0, 0, 0, 0, 0],
     [-1, 1, 0, 0, 0],
     [0, 1/2, -1/2, 0, 0],
     [-1/2, 1, -1/3, -1/6, 0],
     [0, 2/3, -2/3, -1/12, 1/12]]
    >>> res[0]  # no approximation possible, using x_list[0] only
    [0, 0, 0, 0, 0]
    >>> res[1]  # classic forward step approximation
    [-1, 1, 0, 0, 0]
    >>> res[2]  # classic centered approximation
    [0, 1/2, -1/2, 0, 0]
    >>> res[3:]  # higher order approximations
    [[-1/2, 1, -1/3, -1/6, 0], [0, 2/3, -2/3, -1/12, 1/12]]

    Let us compare this to a differently defined x_list. Pay attention to
    foo[i][k] corresponding to the gridpoint defined by x_list[k].

    >>> from sympy import S
    >>> from sympy.calculus import finite_diff_weights
    >>> foo = finite_diff_weights(1, [-S(2), -S(1), S(0), S(1), S(2)], 0)[1]
    >>> foo
    [[0, 0, 0, 0, 0],
     [-1, 1, 0, 0, 0],
     [1/2, -2, 3/2, 0, 0],
     [1/6, -1, 1/2, 1/3, 0],
     [1/12, -2/3, 0, 2/3, -1/12]]
    >>> foo[1]  # not the same and of lower accuracy as res[1]!
    [-1, 1, 0, 0, 0]
    >>> foo[2]  # classic double backward step approximation
    [1/2, -2, 3/2, 0, 0]
    >>> foo[4]  # the same as res[4]
    [1/12, -2/3, 0, 2/3, -1/12]

    Note that, unless you plan on using approximations based on subsets of
    x_list, the order of gridpoints does not matter.


    The capability to generate weights at arbitrary points can be
    used e.g. to minimize Runge's phenomenon by using Chebyshev nodes:

    >>> from sympy import cos, symbols, pi, simplify
    >>> from sympy.calculus import finite_diff_weights
    >>> N, (h, x) = 4, symbols('h x')
    >>> x_list = [x+h*cos(i*pi/(N)) for i in range(N,-1,-1)] # chebyshev nodes
    >>> print(x_list)
    [-h + x, -sqrt(2)*h/2 + x, x, sqrt(2)*h/2 + x, h + x]
    >>> mycoeffs = finite_diff_weights(1, x_list, 0)[1][4]
    >>> [simplify(c) for c in  mycoeffs] #doctest: +NORMALIZE_WHITESPACE
    [(h**3/2 + h**2*x - 3*h*x**2 - 4*x**3)/h**4,
    (-sqrt(2)*h**3 - 4*h**2*x + 3*sqrt(2)*h*x**2 + 8*x**3)/h**4,
    6*x/h**2 - 8*x**3/h**4,
    (sqrt(2)*h**3 - 4*h**2*x - 3*sqrt(2)*h*x**2 + 8*x**3)/h**4,
    (-h**3/2 + h**2*x + 3*h*x**2 - 4*x**3)/h**4]

    Notes
    =====

    If weights for a finite difference approximation of 3rd order
    derivative is wanted, weights for 0th, 1st and 2nd order are
    calculated "for free", so are formulae using subsets of x_list.
    This is something one can take advantage of to save computational cost.
    Be aware that one should define x_list from nearest to farest from
    x_list. If not, subsets of x_list will yield poorer approximations,
    which might not grand an order of accuracy of len(x_list) - order.

    See also
    ========

    sympy.calculus.finite_diff.apply_finite_diff


    References
    ==========

    .. [1] Generation of Finite Difference Formulas on Arbitrarily Spaced
            Grids, Bengt Fornberg; Mathematics of computation; 51; 184;
            (1988); 699-706; doi:10.1090/S0025-5718-1988-0935077-0

    """
    # The notation below closely corresponds to the one used in the paper.
    if order < 0:
        raise ValueError("Negative derivative order illegal.")
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


def apply_finite_diff(order, x_list, y_list, x0=S(0)):
    """
    Calculates the finite difference approximation of
    the derivative of requested order at x0 from points
    provided in x_list and y_list.

    Parameters
    ==========

    order: int
        order of derivative to approximate. 0 corresponds to interpolation.
    x_list: sequence
        Sequence of (unique) values for the independent variable.
    y_list: sequence
        The function value at corresponding values for the independent
        variable in x_list.
    x0: Number or Symbol
        At what value of the independent variable the derivative should be
        evaluated. Defaults to S(0).

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
(-x[i - 1] + x[i])*y[i + 1]/((-x[i - 1] + x[i + 1])*(x[i + 1] - x[i])) - \
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


def as_finite_diff(derivative, points=1, x0=None, wrt=None):
    """
    Returns an approximation of a derivative of a function in
    the form of a finite difference formula. The expression is a
    weighted sum of the function at a number of discrete values of
    (one of) the independent variable(s).

    Parameters
    ==========

    derivative: a Derivative instance (needs to have an variables
        and expr attribute).

    points: sequence or coefficient, optional
        If sequence: discrete values (length >= order+1) of the
        independent variable used for generating the finite
        difference weights.
        If it is a coefficient, it will be used as the step-size
        for generating an equidistant sequence of length order+1
        centered around x0. default: 1 (step-size 1)

    x0: number or Symbol, optional
        the value of the independent variable (wrt) at which the
        derivative is to be approximated. default: same as wrt

    wrt: Symbol, optional
        "with respect to" the variable for which the (partial)
        derivative is to be approximated for. If not provided it
        is required that the Derivative is ordinary. default: None


    Examples
    ========

    >>> from sympy import symbols, Function, exp, sqrt, Symbol, as_finite_diff
    >>> x, h = symbols('x h')
    >>> f = Function('f')
    >>> as_finite_diff(f(x).diff(x))
    -f(x - 1/2) + f(x + 1/2)

    The default step size and number of points are 1 and ``order + 1``
    respectively. We can change the step size by passing a symbol
    as a parameter:

    >>> as_finite_diff(f(x).diff(x), h)
    -f(-h/2 + x)/h + f(h/2 + x)/h

    We can also specify the discretized values to be used in a sequence:

    >>> as_finite_diff(f(x).diff(x), [x, x+h, x+2*h])
    -3*f(x)/(2*h) + 2*f(h + x)/h - f(2*h + x)/(2*h)

    The algorithm is not restricted to use equidistant spacing, nor
    do we need to make the approximation around x0, but we can get
    an expression estimating the derivative at an offset:

    >>> e, sq2 = exp(1), sqrt(2)
    >>> xl = [x-h, x+h, x+e*h]
    >>> as_finite_diff(f(x).diff(x, 1), xl, x+h*sq2)
    2*h*((h + sqrt(2)*h)/(2*h) - (-sqrt(2)*h + h)/(2*h))*f(E*h + x)/\
((-h + E*h)*(h + E*h)) + (-(-sqrt(2)*h + h)/(2*h) - \
(-sqrt(2)*h + E*h)/(2*h))*f(-h + x)/(h + E*h) + \
(-(h + sqrt(2)*h)/(2*h) + (-sqrt(2)*h + E*h)/(2*h))*f(h + x)/(-h + E*h)

    Partial derivatives are also supported:

    >>> y = Symbol('y')
    >>> d2fdxdy=f(x,y).diff(x,y)
    >>> as_finite_diff(d2fdxdy, wrt=x)
    -f(x - 1/2, y) + f(x + 1/2, y)

    See also
    ========

    sympy.calculus.finite_diff.apply_finite_diff
    sympy.calculus.finite_diff.finite_diff_weights

    """
    if wrt is None:
        wrt = derivative.variables[0]
        # we need Derivative to be univariate to guess wrt
        if any(v != wrt for v in derivative.variables):
            raise ValueError('if the function is not univariate' +
                             ' then `wrt` must be given')

    order = derivative.variables.count(wrt)

    if x0 is None:
        x0 = wrt

    if not iterable(points):
        # points is simply the step-size, let's make it a
        # equidistant sequence centered around x0
        if order % 2 == 0:
            # even order => odd number of points, grid point included
            points = [x0 + points*i for i
                      in range(-order//2, order//2 + 1)]
        else:
            # odd order => even number of points, half-way wrt grid point
            points = [x0 + points*i/S(2) for i
                      in range(-order, order + 1, 2)]

    if len(points) < order+1:
        raise ValueError("Too few points for order %d" % order)
    return apply_finite_diff(order, points, [
        derivative.expr.subs({wrt: x}) for x in points], x0)
