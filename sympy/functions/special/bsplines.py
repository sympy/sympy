from sympy.core.basic import Basic, S, C, sympify
from sympy.core.function import expand
from sympy.functions import Piecewise, piecewise_fold
from sympy.functions.elementary.piecewise import ExprCondPair
from sympy.core.sets import Interval


def _add_splines(c, b1, d, b2):
    """Construct c*b1 + d*b2."""
    if b1 == S.Zero or c == S.Zero:
        return expand(piecewise_fold(d*b2))
    if b2 == S.Zero or d == S.Zero:
        return expand(piecewise_fold(c*b1))
    new_args = []
    n_intervals = len(b1.args)
    assert(n_intervals==len(b2.args))
    new_args.append((expand(c*b1.args[0].expr), b1.args[0].cond))
    for i in range(1, n_intervals-1):
        new_args.append((
            expand(c*b1.args[i].expr+d*b2.args[i-1].expr),
            b1.args[i].cond
        ))
    new_args.append((expand(d*b2.args[-2].expr), b2.args[-2].cond))
    new_args.append(b2.args[-1])
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
        Piecewise((1, [0, 1]), (0, True))

    For a given (d, knots) there are len(knots)-d-1 B-splines defined, that
    are indexed by n (starting at 0).

    Here is an example of a cubic B-spline:

        >>> bspline_basis(3, range(5), 0, x)
        Piecewise((x**3/6, [0, 1)), (2/3 - 2*x + 2*x**2 - x**3/2, [1, 2)), (-22/3 + 10*x - 4*x**2 + x**3/2, [2, 3)), (32/3 - 8*x + 2*x**2 - x**3/6, [3, 4]), (0, True))

    By repeating knot points, you can introduce discontinuities in the
    B-splines and their derivatives:

        >>> d = 1
        >>> knots = [0,0,2,3,4]
        >>> bspline_basis(d, knots, 0, x)
        Piecewise((1 - x/2, [0, 2]), (0, True))

    It is quite time consuming to construct and evaluate B-splines. If you
    need to evaluate a B-splines many times, it is best to lambdify them
    first:

        >>> from sympy import lambdify
        >>> d = 3
        >>> knots = range(10)
        >>> b0 = bspline_basis(d, knots, 0, x)
        >>> f = lambdify(x, b0)
        >>> y = f(0.5)

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
            (S.One, Interval(knots[n], knots[n+1], False, True)),
            (0, True)
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
        final_ec_pair = result.args[-2]
        final_cond = final_ec_pair.cond
        final_expr = final_ec_pair.expr
        new_args = final_cond.args[:3] + (False,)
        new_ec_pair = ExprCondPair(final_expr, Interval(*new_args))
        new_args = result.args[:-2] + (new_ec_pair, result.args[-1])
        result = Piecewise(*new_args)
    return result

def bspline_basis_set(d, knots, x):
    """Return the len(knots)-d-1 B-splines at x of degree d with knots.

    This function returns a list of Piecewise polynomials that are the
    len(knots)-d-1 B-splines of degree d for the given knots. This function
    calls bspline_basis(d, knots, n, x) for different values of n.

        >>> from sympy import bspline_basis_set
        >>> from sympy.abc import x
        >>> d = 2
        >>> knots = range(5)
        >>> splines = bspline_basis_set(d, knots, x)
        >>> splines
        [Piecewise((x**2/2, [0, 1)), (-3/2 + 3*x - x**2, [1, 2)), (9/2 - 3*x + x**2/2, [2, 3]), (0, True)), Piecewise((1/2 - x + x**2/2, [1, 2)), (-11/2 + 5*x - x**2, [2, 3)), (8 - 4*x + x**2/2, [3, 4]), (0, True))]
    """
    splines = []
    n_splines = len(knots)-d-1
    return [bspline_basis(d, knots, i, x) for i in range(n_splines)]

