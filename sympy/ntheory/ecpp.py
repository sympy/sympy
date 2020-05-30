#----------------------------------------------------------------------------#
#                                                                            #
#                        Elliptic Curve Primality test                       #
#                                                                            #
#----------------------------------------------------------------------------#



def hilbert(D):
    """
    Returns the class number and the Hilbert class polynomial.

    Parameters
    ==========

    D : Negative fundamental discriminant of the curve

    References
    ==========

    .. [1]  Carl Pomerance and Richard Crandall "Prime Numbers:
        A Computational Perspective" (2nd Ed.), page 360
    """
    from sympy import floor, Symbol, poly, re
    from mpmath import j, qp, exp, pi, fabs, sqrt
    if D >= 0:
        raise ValueError("Only accept negative Discriminant")
    x = Symbol('x')
    T = 1
    b = D % 2
    r = floor(sqrt(abs(D) / 3))
    h = 0
    red = set()
    while(b <= r):
        m = (b*b - D) / 4
        a = 1
        while(a*a <= m):
            if m % a != 0:
                a += 1
                continue
            c = m / a
            if b > a:
                a += 1
                continue
            tau = (-b + j*sqrt(abs(D))) / (2*a)
            tau1 = exp(4*pi*j*tau)
            tau2 = exp(2*pi*j*tau)
            tau1 = qp(tau1)**24 * tau1
            tau2 = qp(tau2)**24 * tau2
            f = tau1 / tau2
            J = (256*f + 1)**3 / f
            if b == a or c == a or b == 0:
                T = T * (x - J)
                h += 1
                red.add((a, b, c))
            else:
                T = T * (x*x - 2*J.real*x + fabs(J)**2)
                h += 2
                red.add((a, b, c))
                red.add((a, -1*b, c))
            a += 1
        b += 2
    new_pol = 0
    coeff = poly(T, x).all_coeffs()
    for i in range(len(coeff) - 1, -1, -1):
        new_pol += int(round(re(coeff[len(coeff) - i - 1]))) * x**i
    return h, poly(new_pol, x), red


def c_smith(p, D):
    """
    Returns solutions (x, y) for the equation
    x**2 + |D|y**2 = 4p

    Parameters
    ==========

    D : Negative integer satisfing -4p < D < 0
    p : odd prime number which is parameter for
        right side of the equation

    References
    ==========

    .. [1]  Carl Pomerance and Richard Crandall "Prime Numbers:
        A Computational Perspective" (2nd Ed.), page 107
    """
    from sympy.ntheory.primetest import is_square
    from sympy.ntheory import jacobi_symbol, sqrt_mod
    from sympy import sqrt, floor

    #Checking bound of D
    if D >= 0 or D <= -4*p:
        raise ValueError("D should be in range -4*p < D < 0")
    if D % 4 not in [0, 1]:
        raise ValueError("D = 0, 1 mod4")

    #If no solution exist
    if p == 2:
        if is_square(D + 8):
            return sqrt(D + 8), 1
        return None
    if jacobi_symbol(D, p) < 1:
        return None

    x_0 = sqrt_mod(D, p)
    if x_0 % 2 != D % 2:
        x_0 = p - x_0
    a, b = 2*p, x_0
    c = floor(2*sqrt(p))

    while b > c:
        a, b = b, a % b

    t = 4*p - b*b
    D = abs(D)
    if t % D != 0:
        return None
    if not is_square(t // D):
        return None
    sols = [(b, sqrt(t // D))]
    sols.append((-sols[0][0], -sols[0][1]))
    return sols
