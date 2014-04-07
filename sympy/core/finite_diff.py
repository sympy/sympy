"""
The method of generating weights for finite differences here
is taken from

    Generation of Finite Difference Formulas on Arbitrarily Spaced Grids, Bengt Fornberg
    Mathematics of compuation, 51, 184, 1988, 699-706 doi:10.1090/S0025-5718-1988-0935077-0

"""

from .singleton import S

def finite_diff_weights(order, x_list, x0):
    """
    Calculates the coefficients for the arbitrarily spaced grid
    provided in x_list for computation of `order`:th derivative
    at x0.

    Algorithm from:
    Generation of Finite Difference Formulas on Arbitrarily Spaced Grids, Bengt Fornberg
    Mathematics of compuation, 51, 184, 1988, 699-706 doi:10.1090/S0025-5718-1988-0935077-0

    The notation below colosely corresponds that used in the paper.
    """
    M = order
    N = len(x_list) - 1
    delta = [[[0 for nu in range(N+1)] for n in range(N+1)] for m in range(M+1)]
    delta[0][0][0]= S(1)
    c1 = S(1)
    for n in range(1,N+1):
        c2 = S(1)
        for nu in range(0,n):
            c3 = x_list[n]-x_list[nu]
            c2 = c2 * c3
            if n <= M: delta[n][n-1][nu]=0
            for m in range(0,min(n,M)+1):
                delta[m][n][nu] = (x_list[n]-x0)*delta[m][n-1][nu] -\
                    m*delta[m-1][n-1][nu]
                delta[m][n][nu] /= c3
        for m in range(0,min(n,M)+1):
            delta[m][n][n] = c1/c2*(m*delta[m-1][n-1][n-1] \
                - (x_list[n-1]-x0)*delta[m][n-1][n-1] )
        c1 = c2
    return delta


def apply_finite_difference(order, x_list, y_list, x0):
    """
    Computes the `order`:th derivative of y
    with respect to x at x0. x_list and y_list
    contains the data points to be used, i.e.
    only supply so many points you think makes sense
    to around x0 when extracting the derivative.

    Order = 0 corresponds to intrpolation.

    Notes
    -----
    Algorithm from:
    Generation of Finite Difference Formulas on Arbitrarily Spaced Grids, Bengt Fornberg
    Mathematics of compuation, 51, 184, 1988, 699-706
    doi:10.1090/S0025-5718-1988-0935077-0

    In the original paper the following holds for the notation
    M = order
    N = len(x_list) - 1
    """
    N = len(x_list) - 1

    if len(x_list) != len(y_list):
        raise ValueError("x_list and y_list not equal in length.")

    delta = finite_diff_weights(order, x_list, x0)

    derivative = 0
    for nu in range(0,len(x_list)):
        derivative += delta[order][N][nu]*y_list[nu]

    return derivative
