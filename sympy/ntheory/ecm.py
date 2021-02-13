from sympy.ntheory import sieve, isprime
from sympy.core.power import integer_log
from sympy.core.compatibility import as_int
import random

rgen = random.Random()

#----------------------------------------------------------------------------#
#                                                                            #
#                   Lenstra's Elliptic Curve Factorization                   #
#                                                                            #
#----------------------------------------------------------------------------#


class Point:
    """Montgomery form of Points in an elliptic curve.
    In this form, the addition and doubling of points
    does not need any y-coordinate information thus
    decreasing the number of operations.
    Using Montgomery form we try to perform point addition
    and doubling in least amount of multiplications.

    The elliptic curve used here is of the form
    (E : b*y**2*z = x**3 + a*x**2*z + x*z**2).
    The a_24 parameter is equal to (a + 2)/4.

    References
    ==========

    .. [1]  http://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
    """

    def __init__(self, x_cord, z_cord, a_24, mod):
        """
        Initial parameters for the Point class.

        Parameters
        ==========

        x_cord : X coordinate of the Point
        z_cord : Z coordinate of the Point
        a_24 : Parameter of the elliptic curve in Montgomery form
        mod : modulus
        """
        self.x_cord = x_cord
        self.z_cord = z_cord
        self.a_24 = a_24
        self.mod = mod

    def __eq__(self, other):
        """Two points are equal if X/Z of both points are equal
        """
        from sympy import mod_inverse
        if self.a_24 != other.a_24 or self.mod != other.mod:
            return False
        return self.x_cord * mod_inverse(self.z_cord, self.mod) % self.mod ==\
            other.x_cord * mod_inverse(other.z_cord, self.mod) % self.mod

    def add(self, Q, diff):
        """
        Add two points self and Q where diff = self - Q. Moreover the assumption
        is self.x_cord*Q.x_cord*(self.x_cord - Q.x_cord) != 0. This algorithm
        requires 6 multiplications. Here the difference between the points
        is already known and using this algorihtm speeds up the addition
        by reducing the number of multiplication required. Also in the
        mont_ladder algorithm is constructed in a way so that the difference
        between intermediate points is always equal to the initial point.
        So, we always know what the difference between the point is.


        Parameters
        ==========

        Q : point on the curve in Montgomery form
        diff : self - Q

        Examples
        ========

        >>> from sympy.ntheory.ecm import Point
        >>> p1 = Point(11, 16, 7, 29)
        >>> p2 = Point(13, 10, 7, 29)
        >>> p3 = p2.add(p1, p1)
        >>> p3.x_cord
        23
        >>> p3.z_cord
        17
        """
        u = (self.x_cord - self.z_cord)*(Q.x_cord + Q.z_cord)
        v = (self.x_cord + self.z_cord)*(Q.x_cord - Q.z_cord)
        add, subt = u + v, u - v
        x_cord = diff.z_cord * add * add % self.mod
        z_cord = diff.x_cord * subt * subt % self.mod
        return Point(x_cord, z_cord, self.a_24, self.mod)

    def double(self):
        """
        Doubles a point in an elliptic curve in Montgomery form.
        This algorithm requires 5 multiplications.

        Examples
        ========

        >>> from sympy.ntheory.ecm import Point
        >>> p1 = Point(11, 16, 7, 29)
        >>> p2 = p1.double()
        >>> p2.x_cord
        13
        >>> p2.z_cord
        10
        """
        u, v = self.x_cord + self.z_cord, self.x_cord - self.z_cord
        u, v = u*u, v*v
        diff = u - v
        x_cord = u*v % self.mod
        z_cord = diff*(v + self.a_24*diff) % self.mod
        return Point(x_cord, z_cord, self.a_24, self.mod)

    def mont_ladder(self, k):
        """
        Scalar multiplication of a point in Montgomery form
        using Montgomery Ladder Algorithm.
        A total of 11 multiplications are required in each step of this
        algorithm.

        Parameters
        ==========

        k : The positive integer multiplier

        Examples
        ========

        >>> from sympy.ntheory.ecm import Point
        >>> p1 = Point(11, 16, 7, 29)
        >>> p3 = p1.mont_ladder(3)
        >>> p3.x_cord
        23
        >>> p3.z_cord
        17
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


def _ecm_one_factor(n, B1=10000, B2=100000, max_curve=200):
    """Returns one factor of n using
    Lenstra's 2 Stage Elliptic curve Factorization
    with Suyama's Parameterization. Here Montgomery
    arithmetic is used for fast computation of addition
    and doubling of points in elliptic curve.

    This ECM method considers elliptic curves in Montgomery
    form (E : b*y**2*z = x**3 + a*x**2*z + x*z**2) and involves
    elliptic curve operations (mod N), where the elements in
    Z are reduced (mod N). Since N is not a prime, E over FF(N)
    is not really an elliptic curve but we can still do point additions
    and doubling as if FF(N) was a field.

    Stage 1 : The basic algorithm involves taking a random point (P) on an
    elliptic curve in FF(N). The compute k*P using Montgomery ladder algorithm.
    Let q be an unknown factor of N. Then the order of the curve E, |E(FF(q))|,
    might be a smooth number that divides k. Then we have k = l * |E(FF(q))|
    for some l. For any point belonging to the curve E, |E(FF(q))|*P = O,
    hence k*P = l*|E(FF(q))|*P. Thus kP.z_cord = 0 (mod q), and the unknownn
    factor of N (q) can be recovered by taking gcd(kP.z_cord, N).

    Stage 2 : This is a continuation of Stage 1 if k*P != O. The idea utilize
    the fact that even if kP != 0, the value of k might miss just one large
    prime divisor of |E(FF(q))|. In this case we only need to compute the
    scalar multiplication by p to get p*k*P = O. Here a second bound B2
    restrict the size of possible values of p.

    Parameters
    ==========

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

    while(curve <= max_curve):
        curve += 1

        #Suyama's Paramatrization
        sigma = rgen.randint(6, n - 1)
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
        g = gcd(Q.z_cord, n)

        #Stage 1 factor
        if g != 1 and g != n:
            return g
        #Stage 1 failure. Q.z = 0, Try another curve
        elif g == n:
            continue

        #Stage 2 - Improved Standard Continuation
        S[1] = Q.double()
        S[2] = S[1].double()
        beta[1] = (S[1].x_cord*S[1].z_cord) % n
        beta[2] = (S[2].x_cord*S[2].z_cord) % n

        for d in range(3, D + 1):
            S[d] = S[d - 1].add(S[1], S[d - 2])
            beta[d] = (S[d].x_cord*S[d].z_cord) % n

        g = 1
        B = B1 - 1
        T = Q.mont_ladder(B - 2*D)
        R = Q.mont_ladder(B)

        for r in  range(B, B2, 2*D):
            alpha = (R.x_cord*R.z_cord) % n
            for q in sieve.primerange(r + 2, r + 2*D + 1):
                delta = (q - r) // 2
                f = (R.x_cord - S[d].x_cord)*(R.z_cord + S[d].z_cord) -\
                alpha + beta[delta]
                g = (g*f) % n
            #Swap
            T, R = R, R.add(S[D], T)
        g = gcd(n, g)

        #Stage 2 Factor found
        if g != 1 and g != n:
            return g

    #ECM failed, Increase the bounds
    raise ValueError("Increase the bounds")


def ecm(n, B1=10000, B2=100000, max_curve=200, seed=1234):
    """Performs factorization using Lenstra's Elliptic curve method.

    This function repeatedly calls `ecm_one_factor` to compute the factors
    of n. First all the small factors are taken out using trial division.
    Then `ecm_one_factor` is used to compute one factor at a time.

    Parameters
    ==========

    n : Number to be Factored
    B1 : Stage 1 Bound
    B2 : Stage 2 Bound
    max_curve : Maximum number of curves generated
    seed : Initialize pseudorandom generator

    Examples
    ========

    >>> from sympy.ntheory import ecm
    >>> ecm(25645121643901801)
    {5394769, 4753701529}
    >>> ecm(9804659461513846513)
    {4641991, 2112166839943}
    """
    _factors = set()
    for prime in sieve.primerange(1, 100000):
        if n % prime == 0:
            _factors.add(prime)
            while(n % prime == 0):
                n //= prime
    rgen.seed(seed)
    while(n > 1):
        try:
            factor = _ecm_one_factor(n, B1, B2, max_curve)
        except ValueError:
            raise ValueError("Increase the bounds")
        _factors.add(factor)
        n //= factor

    factors = set()
    for factor in _factors:
        if isprime(factor):
            factors.add(factor)
            continue
        factors |= ecm(factor)
    return factors
