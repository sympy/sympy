from sympy.ntheory import sieve, isprime
from sympy.core.power import integer_log
from sympy.core.compatibility import as_int
import random

#----------------------------------------------------------------------------#
#                                                                            #
#                   Lenstra's Elliptic Curve Factorization                   #
#                                                                            #
#----------------------------------------------------------------------------#


class Point:
    """Montgomery form of Points in an elliptic curve.
    In this form, the addition and doubling of points
    does not need any y-coordinate information thus
    decreasing the number of operation.

    Parameters:
    ===========

    x_corr : X coordinate of the Point
    z_corr : Z coordinate of the Point
    a_24 : Parameter of the elliptic curve in Montgomery form
    mod : modulus

    References
    ==========

    .. [1]  http://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
    """

    def __init__(self, x_corr, z_corr, a_24, mod):
        self.x_corr = x_corr
        self.z_corr = z_corr
        self.a_24 = a_24
        self.mod = mod

    def add(self, Q, diff):
        """
        Add two points P and Q where diff = P - Q.

        Parameters:
        ===========

        Q : point on the curve in Montgomery form
        diff : P - Q
        """
        u = (self.x_corr - self.z_corr)*(Q.x_corr + Q.z_corr)
        v = (self.x_corr + self.z_corr)*(Q.x_corr - Q.z_corr)
        add, subt = u + v, u - v
        x_corr = diff.z_corr*add*add % self.mod
        z_corr = diff.x_corr*subt*subt % self.mod
        return Point(x_corr, z_corr, self.a_24, self.mod)

    def double(self):
        """
        Doubles a point in an elliptic curve in Montgomery form.
        """
        u, v = self.x_corr + self.z_corr, self.x_corr - self.z_corr
        u, v = u*u, v*v
        diff = u - v
        x_corr = u*v % self.mod
        z_corr = diff*(v + self.a_24*diff) % self.mod
        return Point(x_corr, z_corr, self.a_24, self.mod)

    def mont_ladder(self, k):
        """
        Scalar multiplication of a point in Montgomery form
        using Montgomery Ladder Algorithm.

        Parameters:
        ==========

        k : The positive integer multiplier
        """
        Q = self
        R = self.double()
        for i in bin(k)[3:]:
            if  i  == '1':
                Q = R.add(Q, self)
                R = R.double()
            else:
                R = Q.add(R, self)
                Q = Q.double()
        return Q


def ecm_one_factor(n, B1=10000, B2=100000, max_curve=200, seed=1234):
    """Returns one factor of n using
    Lenstra's 2 stage Elliptic curve Factorization
    with Suyama's Parameterization. Here Montgomery
    arthematic is used for fast computation of addition
    and doubling of points in elliptic curve.

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
    from sympy import gcd, mod_inverse, sqrt
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

    random.seed(seed)
    while(curve <= max_curve):
        curve += 1

        #Suyama's Paramatrization
        sigma = random.randint(6, n)
        u = (sigma*sigma - 5) % n
        v = (4*sigma) % n
        diff = v - u
        u_3 = pow(u, 3, n)

        try:
            C = (pow(diff, 3, n)*(3*u + v)*mod_inverse(4*u_3*v, n) - 2) % n
        except ValueError:
            #If the mod_inverse(4*u_3*v, n) doesn't exist
            return gcd(4*u_3*v, n)

        a24 = (C + 2)*mod_inverse(4, n) % n
        Q = Point(u_3 , pow(v, 3, n), a24, n)
        Q = Q.mont_ladder(k)
        g = gcd(Q.z_corr, n)

        #Stage 1 factor
        if g != 1 and g != n:
            return g
        #Stage 1 failure. Q.z = 0, Try another curve
        elif g == n:
            continue

        #Stage 2 - Improved Standard Continuation
        S[1] = Q.double()
        S[2] = S[1].double()
        beta[1] = (S[1].x_corr*S[1].z_corr) % n
        beta[2] = (S[2].x_corr*S[2].z_corr) % n

        for d in range(3, D + 1):
            S[d] = S[d - 1].add(S[1], S[d - 2])
            beta[d] = (S[d].x_corr*S[d].z_corr) % n

        g = 1
        B = B1 - 1
        T = Q.mont_ladder(B - 2*D)
        R = Q.mont_ladder(B)

        for r in  range(B, B2, 2*D):
            alpha = (R.x_corr*R.z_corr) % n
            for q in sieve.primerange(r + 2, r + 2*D + 1):
                delta = (q - r) // 2
                f = (R.x_corr - S[d].x_corr)*(R.z_corr + S[d].z_corr) - alpha + beta[delta]
                g = (g*f) % n
            #Swap
            TT = R
            R = R.add(S[D], T)
            T = TT
        g = gcd(n, g)

        #Stage 2 Factor found
        if g != 1 and g != n:
            return g

    #ECM failed, Increase the bounds
    raise ValueError("Increase the bounds")


def ecm(n, B1=10000, B2=100000, max_curve=200, increase_bound=False, seed=1234):
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
            a = ecm_one_factor(n, B1, B2, max_curve, seed)
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
