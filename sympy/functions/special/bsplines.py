from sympy.core.basic import Basic, S, C, sympify
from sympy.functions import Piecewise, piecewise_fold
from sympy.core.sets import Interval

def bsplinebasis(x, knots, j, n):
    """The j-th B-spline at x of degree n with knots.

    Everything is indexed starting at zero.
    """
    n_intervals = len(knots)-1
    if n==0:
        args = []
        for i in range(n_intervals):
            if i == j:
                value = S.One
            else:
                value = S.Zero
            args.append( (value, Interval(knots[i], knots[i+1], False, True)) )
        return Piecewise(*args)
    elif n>0:
        result = 0
        denom = knots[j+n] - knots[j]
        if denom != 0:
            coef = (x - knots[j])/denom
            result += coef*bsplinebasis(x, knots, j, n-1)
        denom = knots[j+n+1] - knots[j+1]
        if denom != 0:
            coef = (knots[j+n+1] - x)/denom
            result += coef*bsplinebasis(x, knots, j+1, n-1)
        return piecewise_fold(result)
    else:
        raise ValueError('degree must be non-negative: %r' % n)