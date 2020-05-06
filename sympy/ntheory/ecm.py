class point:
    """For denoting the point on the curve.
    For ECM we don't need the y-coordinate. All
    arithmetic operations use only the x and z coordinates.
    """

    def __init__(self, x, z):
        self.x = x
        self.z = z

def _pt_add(P, Q, diff, n):
    """
    Add two points P and Q where diff = P - Q.
    Parameters:
    ===========
    P : point on the cuvve in Montgomery form
    Q : point on the cuvve in Montgomery form
    diff = P - Q
    n = modulus
    """

    u = (P.x - P.z)*(Q.x + Q.z)
    v = (P.x + P.z)*(Q.x - Q.z)
    s, d = u + v, u - v
    x = diff.z*s*s % n
    z = diff.x*d*d % n
    return point(x, z)

def _pt_double(P, n, a24):
    """
    Doubles a point in an elliptic curve in Montgomery form.
    Parameter:
    P : Point on the curve
    a24 : Parameter for doubling
    n : modulus
    """

    u, v = P.x + P.z, P.x - P.z
    u, v = u*u, v*v
    t = u - v
    x = u*v % n
    z = t*(v + a24*t) % n
    return point(x, z)

def _mont_ladder(k, P, n, a24):
    """
    Scalar multiblication of a point in
    Ellipic Curve in Montgomery form.
    Parameters:
    k : Integer
    P : Poin on the curve in Montgomery form
    n : modulus
    a24 : Parameter for doubling
    """

    Q = P
    R = _pt_double(P, n, a24)
    for i in bin(k)[3:]:
        if  i  == '1':
            Q = _pt_add(R, Q, P, n)
            R = _pt_double(R, n, a24)
        else:
            R = _pt_add(Q, R, P, n)
            Q = _pt_double(Q, n, a24)
    return Q

def ecm_one_factor(n, B1= 10000, B2= 100000, max_curve = 200):
    """Returns one factor of n using
    Lenstra's 2 stage Elliptic curve Factorization
    with Suyama's Parameterization.
    References
    ==========
    .. [1]  Carl Pomerance and Richard Crandall "Prime Numbers:
        A Computational Perspective" (2nd Ed.), page 344
    """

    from sympy.ntheory import sieve, isprime
    from sympy import gcd, mod_inverse, sqrt
    from sympy.core.power import integer_log
    from sympy.core.compatibility import as_int
    import random
    n = as_int(n)
    if B1 % 2 != 0 or B2 % 2 != 0:
        raise ValueError("The Bounds should be an even integer")
    sieve.extend(B2)
    if isprime(n):
        return n
    curve = 0
    D = int(sqrt(B2))
    beta = [0]*(D + 1)
    S = [0]*(D + 1)
    k = 1
    for p in sieve.primerange(1, B1 + 1):
        k *= pow(p, integer_log(B1, p)[0])
    g = 1
    while(curve <= max_curve):
        curve+=1
        #Suyama's Paramatrization
        sigma = random.randint(6,  n )
        u = (sigma*sigma - 5) % n
        v = (4*sigma) % n
        diff = v - u
        u_3 = pow(u, 3, n)
        try:
            A = (pow(diff, 3, n)*(3*u + v)*mod_inverse(4*u_3*v, n) - 2) % n
        except ValueError:
            return g
        a24 = (A + 2)*mod_inverse(4, n) % n
        P = point(u_3 , pow(v, 3, n))
        Q = _mont_ladder(k, P, n, a24)
        g = gcd(n, Q.z)

        #Stage 1 factor
        if g != 1 and g != n:
            return g

        #Stage 2 - Improved Standard Continuation
        S[1] = _pt_double(Q, n, a24)
        S[2] = _pt_double(S[1], n, a24)
        beta[1] = (S[1].x*S[1].z) % n
        beta[2] = (S[2].x*S[2].z) % n
        for d in range(3, D + 1):
            S[d] = _pt_add(S[d - 1], S[1], S[d - 2], n)
            beta[d] = (S[d].x*S[d].z) % n
        g = 1
        B = B1 - 1
        R = _mont_ladder(B, Q, n, a24)
        T = _mont_ladder(B - 2*D, Q, n, a24)
        for r in  range(B, B2, 2*D):
            alpha = (R.x*R.z) % n
            for q in sieve.primerange(r + 2, r + 2*D + 1):
                d = (q - r) // 2
                f = (R.x - S[d].x)*(R.z + S[d].z) - alpha + beta[d]
                g = (g*f) % n
            TT = R
            R = _pt_add(R, S[D], T, n)
            T = TT
        g = gcd(n, g)

        #Stage 2 Factor found
        if g != 1 and g != n:
            return g
    #ECM failed, Increase the bounds
    raise ValueError("Increase the bounds")

def ecm(n):
    """Performs factorization using
    Lenstra's Elliptic curve method
    """
    from sympy import sieve, isprime
    factor = set()
    for i in sieve.primerange(1, 100000):
        if n % i == 0:
            factor.add(i)
            while(n % i == 0):
                n //= i
    while(n % 2 == 0):
        n //= 2
        factor.add(2)
    while(n > 1):
        a = ecm_one_factor(n)
        factor.add(a)
        n //= a
    final = set()
    for i in factor:
        if isprime(i):
            final.add(i)
            continue
        final |= ecm(i)
    return final




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


def cornacchia_smith(p, D):
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
    if t % abs(D) != 0:
        return None
    if not is_square(t // abs(D)):
        return None
    sols = [(b, sqrt(t // abs(D))), (-b, -sqrt(t // abs(D)))]

    return sols
