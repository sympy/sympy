from sympy.core.basic import Basic, S, C, sympify
from sympy.core.function import expand
from sympy.functions import Piecewise, piecewise_fold
from sympy.functions.elementary.piecewise import ExprCondPair
from sympy.core.sets import Interval


def _add_splines(c, b1, d, b2):
    """Construct c*b1 + d*b2."""
    # print 'b1: ', b1
    # print 'b2: ', b2
    if b1 == S.Zero or c == S.Zero:
        return piecewise_fold(expand(d*b2))
    if b2 == S.Zero or d == S.Zero:
        return piecewise_fold(expand(c*b1))
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

    Everything is indexed starting at zero.
    """
    knots = [sympify(k) for k in knots]
    n_knots = len(knots)
    n_intervals = n_knots-1
    if n+d+1 > n_intervals:
        raise ValueError('n+d+1 must not exceed len(knots)-1')
    if d==0:
        result = Piecewise(
            (S.One, Interval(knots[n], knots[n+1], False, True)),
            (0, True)
        )
    elif d>0:
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
    """Return all B-splines at x of degree d with knots.

    >>> d = 2
    >>> knots = range(5)
    >>> splines = bspline_basis_set(d, knots, x)
    >>> splines
    [Piecewise((x**2/2, [0, 1)), (-3/2 + 3*x - x**2, [1, 2)), (9/2 - 3*x + x**2/2, [2, 3]), (0, True)), Piecewise((1/2 - x + x**2/2, [1, 2)), (-11/2 + 5*x - x**2, [2, 3)), (8 - 4*x + x**2/2, [3, 4]), (0, True))]
    """
    splines = []
    n_splines = len(knots)-d-1
    for i in range(n_splines):
        b = bspline_basis(d,knots, i, x)
        splines.append(b)
    return splines
