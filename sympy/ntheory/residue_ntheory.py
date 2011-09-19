from sympy.core.numbers import igcd
from primetest import isprime
from factor_ import factorint, trailing


def totient_(n):
    """returns the number of integers less than n
    and relatively prime to n"""
    if n < 1:
        raise ValueError("n must be a positive integer")
    tot = 0
    for x in xrange(1, n):
        if igcd(x, n) == 1:
            tot += 1
    return tot


def n_order(a, n):
    """ returns the order of a modulo n
    Order of a modulo n is the smallest integer
    k such that a^k leaves a remainder of 1 with n.
    """
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
    >>> set([i**2 % 7 for i in range(7)])
    [0, 1, 2, 4]
    >>> [j for j in range(7) if is_quad_residue(j, 7)]
    [0, 1, 2, 4]
    """
    if any(int(i) != i for i in [a, p]) or p < 1:
        raise ValueError('a and p must be integers and p must be > 1')
    if a == 0 or a == 1 or p == 2:
        return True
    if p == 1:
        return a == 0
    if a > p or a < 0:
        a = a % p
    if not isprime(p):
        if jacobi_symbol(a, p) == -1:
            print 'yo'
            return False
        for i in range(2, p//2):
            if i**2 % p == a:
                return True
        return False

    def square_and_multiply(a, n, p):
        if n == 0:
            return 1
        elif n == 1:
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
    Return 1 if x**2 cong m mod(n) has a solution else return -1.

    jacobi_symbol(m, n) is product of the legendre_symbol(m, p) for all the prime factors p of n.

    **Examples**

    >>> from sympy.ntheory import jacobi_symbol
    >>> jacobi_symbol(45, 77)
    -1
    >>> jacobi_symbol(60, 121)
    1

    """
    if not n % 2:
        raise ValueError("n should be an odd integer")
    if m < 0 or m > n:
        m = m % n
    if igcd(m, n) != 1:
        return 0

    j = 1
    s = trailing(m)
    m = m >> s
    if s % 2 and n % 8 in [3, 5]:
        j = (-1)**s

    while m != 1:
        if m % 4 == 3 and n % 4 == 3:
            j *= -1
        m, n = n % m, m
        s = trailing(m)
        m = m >> s
        if s % 2 and n % 8 in [3, 5]:
            j = (-1)**s
    return j
