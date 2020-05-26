#----------------------------------------------------------------------------#
#                                                                            #
#                   Lenstra's Elliptic Curve Factorization                   #
#                                                                            #
#----------------------------------------------------------------------------#


class Point:
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

    P : point on the curve in Montgomery form
    Q : point on the curve in Montgomery form
    diff = P - Q
    n = modulus
    """

    u = (P.x - P.z)*(Q.x + Q.z)
    v = (P.x + P.z)*(Q.x - Q.z)
    s, d = u + v, u - v
    x = diff.z*s*s % n
    z = diff.x*d*d % n
    return Point(x, z)

def _pt_double(P, n, a24):
    """
    Doubles a point in an elliptic curve in Montgomery form.

    Parameter:
    ==========

    P : Point on the curve
    a24 : Parameter for doubling
    n : modulus
    """

    u, v = P.x + P.z, P.x - P.z
    u, v = u*u, v*v
    t = u - v
    x = u*v % n
    z = t*(v + a24*t) % n
    return Point(x, z)

def _mont_ladder(k, P, n, a24):
    """
    Scalar multiblication of a point in
    Ellipic Curve in Montgomery form.

    Parameters:
    ==========

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

def ecm_one_factor(n, B1=10000, B2=100000, max_curve=200):
    """Returns one factor of n using
    Lenstra's 2 stage Elliptic curve Factorization
    with Suyama's Parameterization.

    Parameters:
    ===========

    n : Number to be Factored
    B1 : Stage 1 Bound
    B2 : Stage 2 Bound
    max_curve : Maximum number of curves generated

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
        curve += 1

        #Suyama's Paramatrization
        sigma = random.randint(6, n)
        u = (sigma*sigma - 5) % n
        v = (4*sigma) % n
        diff = v - u
        u_3 = pow(u, 3, n)

        try:
            A = (pow(diff, 3, n)*(3*u + v)*mod_inverse(4*u_3*v, n) - 2) % n
        except ValueError:
            return g
        a24 = (A + 2)*mod_inverse(4, n) % n
        P = Point(u_3 , pow(v, 3, n))
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

def ecm(n, B1=10000, B2=100000, max_curve=200, increase_bound=False):
    """Performs factorization using
    Lenstra's Elliptic curve method

    Parameters:
    ===========

    n : Number to be Factored
    B1 : Stage 1 Bound
    B2 : Stage 2 Bound
    max_curve : Maximum number of curves generated
    increase_bound : If True, then the Stage 1 and 2 bounds increase
        by a factor of 10 incase of ECM failure
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
        try:
            a = ecm_one_factor(n, B1, B2, max_curve)
        except ValueError:
            if increase_bound:
                B1 *= 10
                B2 *= 10
            else:
                raise ValueError("Increase the bounds")
            continue
        factor.add(a)
        n //= a
    final = set()
    for i in factor:
        if isprime(i):
            final.add(i)
            continue
        final |= ecm(i)
    return final
