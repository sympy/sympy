from sympy.core import Basic
from sympy.core.numbers import igcd
from sympy.utilities.iterables import uniquify

from primetest import isprime
from factor_ import factorint, trailing


def int_tested(*j):
    "Return all args as integers after confirming that they are integers."
    i = tuple([int(i) for i in j])
    if i != j:
        raise ValueError('all arguments were not integers')
    if len(i) == 1:
        return i[0]
    return i

class Residue(Basic):
    """
    The residue classes of a function f(x) mod n are all
    possible values of the residue f(x) (mod n).
    For example, the residue classes of x**2 (mod 6) are
    [0, 1, 3, 4].

    Examples:
    >>> from sympy.ntheory.residue_ntheory import Residue
    >>> from sympy.abc import x
    >>> a = x**2
    >>> b = Residue(a, 6)
    >>> b.values()
    [0, 1, 3, 4]
    """

    _v = None
    _n = None

    @property
    def v(self):
        return self._v

    @property
    def n(self):
        return self._n

    def __new__(cls, *args):
        ret_obj = Basic.__new__(cls, *args)
        val = args[0]
        n = args[1]
        n = int_tested(n)
        ret_obj._v = val
        ret_obj._n = n
        return ret_obj

    def __mul__(self, other):
        """
        Routine for multiplication of residue classes.

        a * b = (a * b) mod n

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> from sympy.abc import x
        >>> a = x
        >>> b = x
        >>> c = Residue(a, 6)
        >>> d = Residue(b, 6)
        >>> (c * d).values()
        [0, 1, 3, 4]
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class")
        if self.n != other.n:
            raise ValueError("Can't multiply two elements from differnt residue classes")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return Residue((self.v * other.v) % self.n, self.n)
        return Residue(self.v * other.v, self.n)

    def __div__(self, other):
        """
        Routine for division of residue classes.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> from sympy.abc import x
        >>> a = x**17
        >>> b = x**15
        >>> c = Residue(a, 6)
        >>> d = Residue(b, 6)
        >>> (c/d).values()
        [0, 1, 3, 4]
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class")
        if self.n != other.n:
            raise ValueError("Can't divide two elements from differnt residue classes")
        return self * other.inv()

    def __add__(self, other):
        """
        Routine for addition of residue classes.

        a + b = (a + b) mod n

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> from sympy.abc import x
        >>> a = x**2 + 5
        >>> b = -5
        >>> c = Residue(a, 6)
        >>> d = Residue(b, 6)
        >>> (c + d).values()
        [0, 1, 3, 4]
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a Residue class")
        if self.n != other.n:
            raise ValueError("Can't add two elements from differnt residue classes")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return Residue((self.v + other.v) % self.n, self.n)
        return Residue(self.v + other.v, self.n)

    def __sub__(self, other):
        """
        Routine for substraction of residue classes.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> from sympy.abc import x
        >>> a = Residue(x**2, 6)
        >>> b = Residue(5, 6)
        >>> (a - b).values()
        [2, 4, 5]
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a Residue class")
        if self.n != other.n:
            raise ValueError("Can't substract two elements from differnt residue classes")
        return self.__add__(-other)

    def __neg__(self):
        """
        Negates the residue class.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> from sympy.abc import x
        >>> b = Residue(x**2 - 5, 6)
        >>> (-b).values()
        [1, 2, 4]
        """
        return Residue(-self.v, self.n)

    def __pow__(self, n):
        """
        Computes the exponent of the residue class.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> from sympy.abc import x
        >>> a = x
        >>> b = Residue(a, 6)
        >>> (b**2).values()
        [0, 1, 3, 4]
        """
        n = int_tested(n)
        n = n % totient_(self.n)
        new = Residue(1, self.n)
        if n == 0:
            return new
        if isinstance(self.v, int):
            new = Residue(pow(self.v, abs(n), self.n), self.n)
        else:
            new = Residue(pow(self.v, abs(n)), self.n)
        if n < 0:
            new = new.inv()
        return new

    def ord(self):
        """
        Exponent of g: power of g > 0 that equals 1

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> a = Residue(4, 7)
        >>> a.ord()
        3
        """
        i = 1
        if isinstance(self.v, int):
            while (self**i).v != 1:
                i += 1
            return i

    def inv(self):
        """
        Computes the inverse of a residue class.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> a = Residue(4, 7)
        >>> a.inv()
        Residue(2, 7)
        """
        return pow(self, totient_(self.n) - 1)

    def __gte__(self, other):
        """
        Checks if a residue class is greater than or equal to another.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> a = Residue(4, 7)
        >>> b = Residue(2, 7)
        >>> a >= b
        False
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class.")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return self.v <= other.v

    def __lte__(self, other):
        """
        Checks if a residue class is greater than or equal to another.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> a = Residue(4, 7)
        >>> b = Residue(2, 7)
        >>> a <= b
        True
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class.")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return self.v >= other.v

    def __lt__(self, other):
        """
        Checks if a residue class is greater than or equal to another.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> a = Residue(4, 7)
        >>> b = Residue(2, 7)
        >>> a < b
        False
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class.")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return self.v < other.v

    def __gt__(self, other):
        """
        Checks if a residue class is greater than or equal to another.

        Examples:
        >>> from sympy.ntheory.residue_ntheory import Residue
        >>> a = Residue(4, 7)
        >>> b = Residue(2, 7)
        >>> a > b
        True
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class.")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return self.v > other.v

    def __eq__(self, other):
        """
        Checks if a residue class is greater than or equal to another.
        """
        if not isinstance(other, Residue):
            raise ValueError("The second operand is not a residue class.")
        if isinstance(self.v, int) and isinstance(other.v, int):
            return self.v == other.v
        else:
            return self.values == other.values

    def values(self):
        """
        List the members of the residue class.
        """
        if isinstance(self._v, int):
            return [self.v % self.n]
        symbol = self._v.as_poly().gens[0]
        ret_obj = sorted(uniquify([self._v.subs(symbol, i) % self.n \
                                   for i in xrange(1, self.n)]))
        if self._v.subs(self._v.args[0], 0) == 0:
            ret_obj = [0] + ret_obj
        return ret_obj

    def __repr__(self):
        return str(self.values())


def totient_(n):
    """returns the number of integers less than n
    and relatively prime to n"""
    n = int_tested(n)
    if n < 1:
        raise ValueError("n must be a positive integer")
    tot = 0
    for x in xrange(1, n):
        if igcd(x, n) == 1:
            tot += 1
    return tot

def relprimes(n):
    """
    List integers 0 < i < n, such that i, n are coprime
    """
    return [x for x in xrange(1, n) if igcd(n,x) == 1]

def n_order(a, n):
    """ returns the order of a modulo n
    Order of a modulo n is the smallest integer
    k such that a^k leaves a remainder of 1 with n.
    """
    a, n = int_tested(a, n)
    if igcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    group_order = totient_(n)
    factors = factorint(group_order)
    order = 1
    if a > n:
        a = a % n
    for p, e in factors.iteritems():
        exponent = group_order
        for f in xrange(0, e + 1):
            if (a ** (exponent)) % n != 1:
                order *= p ** (e - f + 1)
                break
            exponent = exponent // p
    return order


def is_primitive_root(a, p):
    """
    returns True if a is a primitive root of p
    """
    a, p = int_tested(a, p)
    if igcd(a, p) != 1:
        raise ValueError("The two numbers should be relatively prime")
    if a > p:
        a = a % p
    if n_order(a, p) == totient_(p):
        return True
    else:
        return False


def is_quad_residue(a, p):
    """
    Returns True if ``a`` (mod ``p``) is in the set of squares mod ``p``,
    i.e a % p in set([i**2 % p for i in range(p)]). If ``p`` is an odd
    prime, an iterative method is used to make the determination.

    >>> from sympy.ntheory import is_quad_residue
    >>> list(set([i**2 % 7 for i in range(7)]))
    [0, 1, 2, 4]
    >>> [j for j in range(7) if is_quad_residue(j, 7)]
    [0, 1, 2, 4]
    """
    a, p = int_tested(a, p)
    if p < 1:
        raise ValueError('p must be > 0')
    if a >= p or a < 0:
        a = a % p
    if a < 2 or p < 3:
        return True
    if not isprime(p):
        if p % 2 and jacobi_symbol(a, p) == -1:
            return False
        for i in range(2, p//2 + 1):
            if i**2 % p == a:
                return True
        return False

    def square_and_multiply(a, n, p):
        if n == 1:
            return a
        elif n % 2 == 1:
            return ((square_and_multiply(a, n // 2, p) ** 2) * a) % p
        else:
            return (square_and_multiply(a, n // 2, p) ** 2) % p

    return (square_and_multiply(a, (p - 1) // 2, p) % p) == 1


def legendre_symbol(a, p):
    """
    return 1 if a is a quadratic residue of p, 0 if a is multiple of p,
    else return -1
    p should be an odd prime by definition
    """
    a, p = int_tested(a, p)
    if not isprime(p) or p == 2:
        raise ValueError("p should be an odd prime")
    _, a = divmod(a, p)
    if not a:
        return 0
    if is_quad_residue(a, p):
        return 1
    else:
        return -1


def jacobi_symbol(m, n):
    """
    Returns 0 if m cong 0 mod(n),
            1 if x**2 cong m mod(n) has a solution, else
           -1.

    jacobi_symbol(m, n) is product of the legendre_symbol(m, p)
    for all the prime factors p of n.

    **Examples**

    >>> from sympy.ntheory import jacobi_symbol, legendre_symbol
    >>> from sympy import Mul, S
    >>> jacobi_symbol(45, 77)
    -1
    >>> jacobi_symbol(60, 121)
    1

    The relationship between the jacobi_symbol and legendre_symbol can
    be demonstrated as follows:
        >>> L = legendre_symbol
        >>> S(45).factors()
        {3: 2, 5: 1}
        >>> jacobi_symbol(7, 45) == L(7, 3)**2 * L(7, 5)**1
        True
    """
    m, n = int_tested(m, n)
    if not n % 2:
        raise ValueError("n should be an odd integer")
    if m < 0 or m > n:
        m = m % n
    if not m:
        return int(n == 1)
    if n == 1 or m == 1:
        return 1
    if igcd(m, n) != 1:
        return 0

    j = 1
    s = trailing(m)
    m = m >> s
    if s % 2 and n % 8 in [3, 5]:
        j *= -1

    while m != 1:
        if m % 4 == 3 and n % 4 == 3:
            j *= -1
        m, n = n % m, m
        s = trailing(m)
        m = m >> s
        if s % 2 and n % 8 in [3, 5]:
            j *= -1
    return j

def bin_gcd(a, b):
    """
    Extended version of Euclid's Algorithm (binary GCD)
    Returns (m, n, gcd) such that  (m * a) + (n * b) = gcd(a, b)

    Examples:
    >>> from sympy.ntheory.residue_ntheory import bin_gcd
    >>> bin_gcd(10, 15)
    (14, -9, 5)
    """
    g, u, v = [b, a], [1, 0], [0, 1]
    while g[1] != 0:
        y = g[0] // g[1]
        g[0], g[1] = g[1], g[0] % g[1]
        u[0], u[1] = u[1], u[0] - (y * u[1])
        v[0], v[1] = v[1], v[0] - (y * v[1])
    m = v[0] % b
    gcd = (m * a) % b
    n = (gcd - m * a) // b
    return (m, n, gcd)

