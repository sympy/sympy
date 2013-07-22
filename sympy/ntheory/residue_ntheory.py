from sympy.core.numbers import igcd
from sympy.core.compatibility import as_int
from .primetest import isprime
from .factor_ import factorint, trailing, totient


def n_order(a, n):
    """Returns the order of ``a`` modulo ``n``.

    The order of ``a`` modulo ``n`` is the smallest integer
    ``k`` such that ``a**k`` leaves a remainder of 1 with ``n``.

    Examples
    ========

    >>> from sympy.ntheory import n_order
    >>> n_order(3, 7)
    6
    >>> n_order(4, 7)
    3
    """
    a, n = as_int(a), as_int(n)
    if igcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    group_order = totient(n)
    factors = factorint(group_order)
    order = 1
    if a > n:
        a = a % n
    for p, e in factors.items():
        exponent = group_order
        for f in range(e + 1):
            if a**exponent % n != 1:
                order *= p ** (e - f + 1)
                break
            exponent = exponent // p
    return order


def is_primitive_root(a, p):
    """
    Returns True if ``a`` is a primitive root of ``p``

    ``a`` is said to be the primitive root of ``p`` if gcd(a, p) == 1 and
    totient(p) is the smallest positive number s.t.

        a**totient(p) cong 1 mod(p)

    Examples
    ========

    >>> from sympy.ntheory import is_primitive_root, n_order, totient
    >>> is_primitive_root(3, 10)
    True
    >>> is_primitive_root(9, 10)
    False
    >>> n_order(3, 10) == totient(10)
    True
    >>> n_order(9, 10) == totient(10)
    False

    """
    a, p = as_int(a), as_int(p)
    if igcd(a, p) != 1:
        raise ValueError("The two numbers should be relatively prime")
    if a > p:
        a = a % p
    return n_order(a, p) == totient(p)


def is_quad_residue(a, p):
    """
    Returns True if ``a`` (mod ``p``) is in the set of squares mod ``p``,
    i.e a % p in set([i**2 % p for i in range(p)]). If ``p`` is an odd
    prime, an iterative method is used to make the determination:

    >>> from sympy.ntheory import is_quad_residue
    >>> list(set([i**2 % 7 for i in range(7)]))
    [0, 1, 2, 4]
    >>> [j for j in range(7) if is_quad_residue(j, 7)]
    [0, 1, 2, 4]

    See Also
    ========

    legendre_symbol, jacobi_symbol
    """
    a, p = as_int(a), as_int(p)
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

    return pow(a, (p - 1) // 2, p) == 1


def legendre_symbol(a, p):
    """
    Returns
    =======

    1. 0 if a is multiple of p
    2. 1 if a is a quadratic residue of p
    3. -1 otherwise

    p should be an odd prime by definition

    Examples
    ========

    >>> from sympy.ntheory import legendre_symbol
    >>> [legendre_symbol(i, 7) for i in range(7)]
    [0, 1, 1, -1, 1, -1, -1]
    >>> list(set([i**2 % 7 for i in range(7)]))
    [0, 1, 2, 4]

    See Also
    ========

    is_quad_residue, jacobi_symbol

    """
    a, p = as_int(a), as_int(p)
    if not isprime(p) or p == 2:
        raise ValueError("p should be an odd prime")
    a = a % p
    if not a:
        return 0
    if is_quad_residue(a, p):
        return 1
    return -1


def jacobi_symbol(m, n):
    """
    Returns the product of the legendre_symbol(m, p)
    for all the prime factors, p, of n.

    Returns
    =======

    1. 0 if m cong 0 mod(n)
    2. 1 if x**2 cong m mod(n) has a solution
    3. -1 otherwise

    Examples
    ========

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

    See Also
    ========

    is_quad_residue, legendre_symbol
    """
    m, n = as_int(m), as_int(n)
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
