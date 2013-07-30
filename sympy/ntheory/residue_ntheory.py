from __future__ import print_function, division

from sympy.core.numbers import igcd, igcdex
from sympy.core.compatibility import as_int, xrange
from .primetest import isprime
from .factor_ import factorint, trailing, totient
from random import randint
from itertools import product

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
    for p, e in factors.iteritems():
        exponent = group_order
        for f in xrange(e + 1):
            if pow(a, exponent, n) != 1:
                order *= p ** (e - f + 1)
                break
            exponent = exponent // p
    return order

def primitive_root(p, all_roots=False):
    """
    compute the smallest primitive root

    References
    ==========

    [1] W. Stein "Elementary Number Theory" (2011), page 44
    [2] P. Hackman "Elementary Number Theory" (2009),  Chapter C

    Parameters
    ==========

    p : integer
    all_roots : if True the list of all primitive roots is returned

    TODO : case in which ``p`` is not prime

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import primitive_root
    >>> primitive_root(19)
    2
    >>> primitive_root(19, True)
    [2, 3, 10, 13, 14, 15]
    """
    p = as_int(p)
    f = factorint(p)
    if len(f) > 2:
        return None
    if len(f) == 2:
        if 2 not in f or f[2] > 1:
            return None

        # case p = 2*p1**k, p1 prime
        for p1, e1 in f.items():
            if p1 != 2:
                break
        # see Ref [2], page 72
        return primitive_root(p1, all_roots)
    else:
        if 2 in f:
            if p == 2:
                return 1
            if p == 4:
                return 3
            return None
        p1, n = f.items()[0]
        if n > 1:
            # see Ref [2], page 81
            gv = primitive_root(p1, all_roots)
            if not all_roots:
                gv = [gv]
            res = []
            for g in gv:
                if is_primitive_root(g, p1**2):
                    res.append(p)
                else:
                    res.append(g + p1)
            if all_roots:
                res.sort()
                return res
            else:
                return min(res)

    # see [1]
    v = [(p - 1) // i for i in factorint(p - 1).keys()]
    pv = []
    for a in xrange(2, p):
        for pw in v:
            if pow(a, pw, p) == 1:
                break
        else:
            if not all_roots:
                return a
            pv.append(a)
    return pv

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

def _sqrt_mod_tonelli_shanks(a, p):
    """
    Returns the square root in the case of ``p`` prime with ``p == 1 (mod 8)``

    References
    ==========

    R. Crandall and C. Pomerance "Prime Numbers", 2nt Ed., page 101
    """
    s = trailing(p - 1)
    t = p >> s
    # find a non-quadratic residue
    while 1:
        d = randint(2, p - 1)
        r = legendre_symbol(d, p)
        if r == -1:
            break
    #assert legendre_symbol(d, p) == -1
    A = pow(a, t, p)
    D = pow(d, t, p)
    m = 0
    for i in xrange(s):
        adm = A*pow(D, m, p) % p
        adm = pow(adm, 2**(s - 1 - i), p)
        if adm % p == p - 1:
            m += 2**i
    #assert A*pow(D, m, p) % p == 1
    x = pow(a, (t + 1)//2, p)*pow(D, m//2, p) % p
    return x

def sqrt_mod(a, p, all_roots=False):
    """
    find the solutions to ``x**2 = a mod p``

    Parameters
    ==========

    a : integer
    p : positive integer
    all_roots : if False returns the smallest root, else the list of roots

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import sqrt_mod
    >>> sqrt_mod(11, 43)
    21
    >>> sqrt_mod(11, 43, True)
    [21, 22]
    """
    from sympy.ntheory.modular import crt
    a, p = as_int(a), as_int(p)

    if isprime(p):
        res = _sqrt_mod_prime_power(a, p, 1)
    else:
        f = factorint(p)
        v = []
        pv = []
        for px, ex in f.items():
            if a % px == 0:
                rx = _sqrt_mod1(a, px, ex)
            else:
                rx = _sqrt_mod_prime_power(a, px, ex)
            if rx is None:
                return None
            v.append(rx)
            pv.append(px**ex)
        res = []
        for vx in product(*v):
            r = crt(pv, vx)[0]
            if r not in res:
                res.append(r)
            if p - r not in res:
                if r:
                    res.append(p - r)
    if res is None:
        return None
    if all_roots:
        res.sort()
    else:
        res = min(res)
    return res

def _sqrt_mod_prime_power(a, p, k):
    """
    find the solutions to ``x**2 = a mod p**k``

    Parameters
    ==========

    a : integer
    p : prime number
    k : semipositive integer

    Notes
    =====

    The modulus ``m`` is ``p**k`` if ``p`` is prime and ``k > 0``

    References
    ==========

    P. Hackman "Elementary Number Theory" (2009),  page 160

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _sqrt_mod_prime_power
    >>> _sqrt_mod_prime_power(11, 43, 1)
    [21, 22]
    """
    from sympy.core.numbers import igcdex

    pk = p**k
    a = a % pk
    if a == 0:
        res = _sqrt_mod1(a, p, k)
        return res

    if k == 1:
        if p == 2:
            return [a]
        if not is_quad_residue(a, p):
            return None

        if p % 4 == 3:
            res = pow(a, (p + 1) // 4, p)
        elif p % 8 == 5:
            sign = pow(a, (p - 1) // 4, p)
            if sign == 1:
                res = pow(a, (p + 3) // 8, p)
            else:
                b = pow(4*a, (p - 5) // 8, p)
                x =  (2*a*b) % p
                if pow(x, 2, p) == a:
                    res = x
                else:
                    res = None
        else:
            res = _sqrt_mod_tonelli_shanks(a, p)

        if res is not None:
            return [res, p - res]
        return None

    if k > 1:
        f = factorint(a)
        if p in f:
            if f[p] % 2 == 1:
                return None
        if p == 2:
            if a % 8 != 1:
                return None
            if k <= 3:
               s = set()
               for i in xrange(0, pk, 4):
                    s.add(1 + i)
                    s.add(-1 + i)
               return list(s)
            rv = [1, 3, 5, 7]
            n = 3
            res = []
            for r in rv:
                nx = n
                while nx < k:
                    r1 = (r**2 - a) >> nx
                    i = r1 % 2
                    r = r + (i << (nx - 1))
                    nx += 1
                if r not in res:
                    res.append(r)
                if pk - r not in res:
                    res.append(pk - r)
            return res
        rv = _sqrt_mod_prime_power(a, p, 1)
        if not rv:
            return None
        r = rv[0]
        fr = r**2 - a
        # hensel lifting
        n = 1
        px = p
        while 1:
            n1 = n
            n1 *= 2
            if n1 > k:
                break
            n = n1
            px = px**2
            frinv = igcdex(2*r, px)[0]
            r = (r - fr*frinv) % px
            fr = r**2 - a
        if n < k:
            px = p**k
            frinv = igcdex(2*r, px)[0]
            r = (r - fr*frinv) % px
        return [r, px - r]

def _sqrt_mod1(a, p, n):
    """
    case in which ``a % p == 0 and n > 0``

    see http://www.numbertheory.org/php/squareroot.html
    """
    pn = p**n
    if a % pn == 0:
        # case gcd(a, p**k) = p**r, r < n
        m = n // 2
        if n % 2 == 1:
            return xrange(0, pn, p**(m + 1))
        else:
            return xrange(0, pn, p**m)
    # case gcd(a, p**k) = p**n
    f = factorint(a)
    r = f[p]
    if r % 2 == 1:
        return None
    m = r // 2
    a1 = a >> r
    if p == 2:
        if n - r == 1:
            res = xrange(1, 2**(n - m + 1), 2)
            res = [x << m for x in res]
            return res
        if n - r == 2:
            res = _sqrt_mod_prime_power(a1, p, n - r)
            if res is None:
                return None
            res = [x << m for x in res]
            s = set()
            for r in res:
                for i in xrange(0, pn, 2**(n - m)):
                    s.add(r + i)
            return list(s)
        if n - r > 2:
            res = _sqrt_mod_prime_power(a1, p, n - r)
            if res is None:
                return None
            res = [x << m for x in res]
            s = set()
            for x in res:
                for i in xrange(0, pn, 2**(n - m - 1)):
                    s.add((x + i) % pn)
            return list(s)
    else:
        if r % 2 == 1:
            return None
        m = r // 2
        a1 = a // p**r
        res1 = _sqrt_mod_prime_power(a1, p, n - r)
        if res1 is None:
            return None
        s = set()
        for x in res1:
            for i in xrange(0, pn, p**(n-r)):
                s.add((x + i) % pn)
        res = list(s)
        res = [x*p**m for x in res]
        return res


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
            if pow(i, 2, p) == a:
                return True
        return False

    return pow(a, (p - 1) // 2, p) == 1


def is_nthpow_residue(a, n, m):
    """
    Returns True if ``x**n == a (mod m)`` has solutions.

    References
    ==========

    P. Hackman "Elementary Number Theory" (2009),  page 76
    """
    if n == 1:
        return True
    if n == 2:
        return is_quad_residue(a, m)
    f = totient(m)
    k = f // igcd(f, n)
    return pow(a, k, m) == 1

def _nthroot_mod1(s, q, p, all_roots):
    """
    Root of ``x**n = s mod p``, ``p`` prime and ``n`` divides ``p - 1``

    References
    ==========

    [1] A. M. Johnston "A Generalized qth Root Algorithm"
    """
    if not isprime(q):
        raise NotImplementedError
    g = primitive_root(p)
    f = p - 1
    assert (p - 1) % q == 0
    # determine k
    k = 0
    while f % q == 0:
        k += 1
        f = f // q
    # find z, x, r1
    f1 = igcdex(-f, q)[0] % q
    z = f*f1
    x = (1 + z) // q
    w = pow(g, z, p)
    r1 = pow(s, x, p)
    s1 = pow(s, f, p)
    y = pow(g, f, p)
    h = pow(g, f*q, p)
    # find t discrete log of s1 base h, h**x = s1 mod p
    # used a naive implementation
    # TODO implement using Ref [1]
    pr = 1
    for t in xrange(p):
        if pr == s1:
            break
        pr = pr*h % p

    g2 = pow(g, z*t, p)
    g3 = igcdex(g2, p)[0]
    r = r1*g3 % p
    #assert pow(r, q, p) == s
    res = [r]
    h = pow(g, (p - 1) // q, p)
    #assert pow(h, q, p) == 1
    hx = r
    for i in range(q - 1):
        hx = (hx*h) % p
        res.append(hx)
    if all_roots:
        res.sort()
        return res
    return min(res)

def nthroot_mod(a, n, p, all_roots=False):
    """
    find the solutions to ``x**n = a mod p``

    Parameters
    ==========

    a : integer
    n : positive integer
    p : positive integer
    all_roots : if False returns the smallest root, else the list of roots

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import nthroot_mod
    >>> nthroot_mod(11, 4, 19)
    8
    >>> nthroot_mod(11, 4, 19, True)
    [8, 11]
    >>> nthroot_mod(68, 3, 109)
    23
    """
    from sympy.core.numbers import igcdex
    if n == 2:
        return sqrt_mod(a, p, all_roots)
    f = totient(p)
    # see Hackman "Elementary Number Theory" (2009), page 76
    if pow(a, f // igcd(f, n), p) != 1:
        return None
    if not isprime(p):
        raise NotImplementedError

    if (p - 1) % n == 0:
        return _nthroot_mod1(a, n, p, all_roots)
    # The roots of ``x**n - a = 0 (mod p)`` are roots of
    # ``gcd(x**n - a, x**(p - 1) - 1) = 0 (mod p)``
    pa = n
    pb = p - 1
    b = 1
    if pa < pb:
        a, pa, b, pb = b, pb, a, pa
    while pb:
        # x**pa - a = 0; x**pb - b = 0
        # x**pa - a = x**(q*pb + r) - a = (x**pb)**q * x**r - a =
        #             b**q * x**r - a; x**r - c = 0; c = b**-q * a mod p
        q, r = divmod(pa, pb)
        c = pow(b, q, p)
        c = igcdex(c, p)[0]
        c = (c * a) % p
        pa, pb = pb, r
        a, b = b, c
    if pa == 1:
        if all_roots:
            res = [a]
        else:
            res = a
    elif pa == 2:
        res = sqrt_mod(a, p, all_roots)
    else:
        res = _nthroot_mod1(a, pa, p, all_roots)
    return res

def quadratic_residues(p):
    """
    Returns the list of quadratic residues.

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import quadratic_residues
    >>> quadratic_residues(7)
    [0, 1, 2, 4]
    """
    r = set()
    for i in xrange((p + 1) // 2):
        r.add(pow(i, 2, p))
    return sorted(list(r))


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
