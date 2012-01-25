from sympy.core import S, sympify, expand
from sympy.functions import Piecewise, piecewise_fold
from sympy.core.sets import Interval


def _add_splines(c, b1, d, b2):
    """Construct c*b1 + d*b2."""
    if b1 == S.Zero or c == S.Zero:
        return expand(piecewise_fold(d*b2))
    if b2 == S.Zero or d == S.Zero:
        return expand(piecewise_fold(c*b1))
    new_args = []
    n_intervals = len(b1.exprcondpairs)
    assert(n_intervals==len(b2.exprcondpairs))
    new_args.append((expand(c*b1.exprcondpairs[0][0]), b1.exprcondpairs[0][1]))
    for i in range(1, n_intervals):
        new_args.append((
            expand(c*b1.exprcondpairs[i][0]+d*b2.exprcondpairs[i-1][0]),
            b1.exprcondpairs[i][1]
        ))
    new_args.append((expand(d*b2.exprcondpairs[-1][0]), b2.exprcondpairs[-1][1]))
    new_args.append(d*b2.otherwise)
    return Piecewise(*new_args)

def bspline_basis(d, knots, n, x, close=True):
    """The n-th B-spline at x of degree d with knots.

    B-Splines are piecewise polynomials of degree d [1].  They are defined on
    a set of knots, which is a sequence of integers or floats.

    The 0th degree splines have a value of one on a single interval:

        >>> from sympy import bspline_basis
        >>> from sympy.abc import x
        >>> d = 0
        >>> knots = range(5)
        >>> bspline_basis(d, knots, 0, x)
        Piecewise((1, And(0 <= x, x <= 1)), 0)

    For a given (d, knots) there are len(knots)-d-1 B-splines defined, that
    are indexed by n (starting at 0).

    Here is an example of a cubic B-spline:

        >>> bspline_basis(3, range(5), 0, x)
        Piecewise((x**3/6, And(0 <= x, x < 1)), (-x**3/2 + 2*x**2 - 2*x + 2/3, And(1 <= x, x < 2)), (x**3/2 - 4*x**2 + 10*x - 22/3, And(2 <= x, x < 3)), (-x**3/6 + 2*x**2 - 8*x + 32/3, And(3 <= x, x <= 4)), 0)

    By repeating knot points, you can introduce discontinuities in the
    B-splines and their derivatives:

        >>> d = 1
        >>> knots = [0,0,2,3,4]
        >>> bspline_basis(d, knots, 0, x)
        Piecewise((-x/2 + 1, And(0 <= x, x <= 2)), 0)

    It is quite time consuming to construct and evaluate B-splines. If you
    need to evaluate a B-splines many times, it is best to lambdify them
    first:

        >>> from sympy import lambdify
        >>> d = 3
        >>> knots = range(10)
        >>> b0 = bspline_basis(d, knots, 0, x)
        >>> f = lambdify(x, b0)
        >>> y = f(0.5)

    See Also
    ========

    bsplines_basis_set

    References
    ==========

    [1] http://en.wikipedia.org/wiki/B-spline

    """
    knots = [sympify(k) for k in knots]
    d = int(d)
    n = int(n)
    n_knots = len(knots)
    n_intervals = n_knots-1
    if n+d+1 > n_intervals:
        raise ValueError('n+d+1 must not exceed len(knots)-1')
    if d==0:
        result = Piecewise(
            (S.One, Interval(knots[n], knots[n+1], False, True).as_relational(x)),
            0
        )
    elif d > 0:
        denom = knots[n+d] - knots[n]
        if denom != S.Zero:
            A = (x - knots[n])/denom
            b1 = bspline_basis(d-1, knots, n, x, close=False)
        else:
            b1 = A = S.Zero

        denom = knots[n+d+1] - knots[n+1]
        if denom != S.Zero:
            B = (knots[n+d+1] - x)/denom
            b2 = bspline_basis(d-1, knots, n+1, x, close=False)
        else:
            b2 = B = S.Zero
        result = _add_splines(A, b1, B, b2)
    else:
        raise ValueError('degree must be non-negative: %r' % n)
    if close:
        final_ec_pair = result.exprcondpairs[-1]
        final_cond = final_ec_pair[1]
        final_expr = final_ec_pair[0]
        final_cond = final_cond.args
        from sympy import And
        if final_cond[0].lts == x:
            final_cond = And(final_cond[0].lts <= final_cond[0].gts, final_cond[1])
        else:
            final_cond = And(final_cond[0], final_cond[1].lts <= final_cond[1].gts)
        new_ec_pair = (final_expr, final_cond)
        new_args = tuple(result.exprcondpairs[:-1]) + (new_ec_pair, result.otherwise)
        result = Piecewise(*new_args)
    return result

def bspline_basis_set(d, knots, x):
    """Return the len(knots)-d-1 B-splines at x of degree d with knots.

    This function returns a list of Piecewise polynomials that are the
    len(knots)-d-1 B-splines of degree d for the given knots. This function
    calls bspline_basis(d, knots, n, x) for different values of n.

    Examples
    ========

    >>> from sympy import bspline_basis_set
    >>> from sympy.abc import x
    >>> d = 2
    >>> knots = range(5)
    >>> splines = bspline_basis_set(d, knots, x)
    >>> splines
    [Piecewise((x**2/2, And(0 <= x, x < 1)), (-x**2 + 3*x - 3/2, And(1 <= x, x < 2)), (x**2/2 - 3*x + 9/2, And(2 <= x, x <= 3)), 0), Piecewise((x**2/2 - x + 1/2, And(1 <= x, x < 2)), (-x**2 + 5*x - 11/2, And(2 <= x, x < 3)), (x**2/2 - 4*x + 8, And(3 <= x, x <= 4)), 0)]

    See Also
    ========

    bsplines_basis
    """
    n_splines = len(knots)-d-1
    return [bspline_basis(d, knots, i, x) for i in range(n_splines)]

