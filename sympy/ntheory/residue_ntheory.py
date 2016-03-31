# -*- coding: utf-8 -*-

from __future__ import print_function, division

from sympy.core.singleton import S
from sympy.core.numbers import igcd, igcdex, mod_inverse
from sympy.core.compatibility import as_int, range
from sympy.core.function import Function
from .primetest import isprime
from .factor_ import factorint, trailing, totient, multiplicity
from random import randint



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
    from collections import defaultdict
    a, n = as_int(a), as_int(n)
    if igcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    factors = defaultdict(int)
    f = factorint(n)
    for px, kx in f.items():
        if kx > 1:
            factors[px] += kx - 1
        fpx = factorint(px - 1)
        for py, ky in fpx.items():
            factors[py] += ky
    group_order = 1
    for px, kx in factors.items():
        group_order *= px**kx
    order = 1
    if a > n:
        a = a % n
    for p, e in factors.items():
        exponent = group_order
        for f in range(e + 1):
            if pow(a, exponent, n) != 1:
                order *= p ** (e - f + 1)
                break
            exponent = exponent // p
    return order


def _primitive_root_prime_iter(p):
    """
    Generates the primitive roots for a prime ``p``

    References
    ==========

    [1] W. Stein "Elementary Number Theory" (2011), page 44

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _primitive_root_prime_iter
    >>> list(_primitive_root_prime_iter(19))
    [2, 3, 10, 13, 14, 15]
    """
    p = as_int(p)
    v = [(p - 1) // i for i in factorint(p - 1).keys()]
    a = 2
    while a < p:
        for pw in v:
            if pow(a, pw, p) == 1:
                break
        else:
            yield a
        a += 1


def primitive_root(p):
    """
    Returns the smallest primitive root or None

    References
    ==========

    [1] W. Stein "Elementary Number Theory" (2011), page 44
    [2] P. Hackman "Elementary Number Theory" (2009),  Chapter C

    Parameters
    ==========

    p : positive integer

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import primitive_root
    >>> primitive_root(19)
    2
    """
    p = as_int(p)
    if p < 1:
        raise ValueError('p is required to be positive')
    if p <= 2:
        return 1
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
        i = 1
        while i < p:
            i += 2
            if i % p1 == 0:
                continue
            if is_primitive_root(i, p):
                return i

    else:
        if 2 in f:
            if p == 4:
                return 3
            return None
        p1, n = list(f.items())[0]
        if n > 1:
            # see Ref [2], page 81
            g = primitive_root(p1)
            if is_primitive_root(g, p1**2):
                return g
            else:
                for i in range(2, g + p1 + 1):
                    if igcd(i, p) == 1 and is_primitive_root(i, p):
                        return i

    return next(_primitive_root_prime_iter(p))


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
    for i in range(s):
        adm = A*pow(D, m, p) % p
        adm = pow(adm, 2**(s - 1 - i), p)
        if adm % p == p - 1:
            m += 2**i
    #assert A*pow(D, m, p) % p == 1
    x = pow(a, (t + 1)//2, p)*pow(D, m//2, p) % p
    return x


def sqrt_mod(a, p, all_roots=False):
    """
    find a root of ``x**2 = a mod p``

    Parameters
    ==========

    a : integer
    p : positive integer
    all_roots : if True the list of roots is returned or None

    Notes
    =====

    If there is no root it is returned None; else the returned root
    is less or equal to ``p // 2``; in general is not the smallest one.
    It is returned ``p // 2`` only if it is the only root.

    Use ``all_roots`` only when it is expected that all the roots fit
    in memory; otherwise use ``sqrt_mod_iter``.

    Examples
    ========

    >>> from sympy.ntheory import sqrt_mod
    >>> sqrt_mod(11, 43)
    21
    >>> sqrt_mod(17, 32, True)
    [7, 9, 23, 25]
    """
    if all_roots:
        return sorted(list(sqrt_mod_iter(a, p)))
    try:
        p = abs(as_int(p))
        it = sqrt_mod_iter(a, p)
        r = next(it)
        if r > p // 2:
            return p - r
        elif r < p // 2:
            return r
        else:
            try:
                r = next(it)
                if r > p // 2:
                    return p - r
            except StopIteration:
                pass
            return r
    except StopIteration:
        return None


def _product(*iters):
    """
    cartesian product generator

    Notes
    =====

    Unlike itertools.product, it works also with iterables which do not fit
    in memory. See http://bugs.python.org/issue10109

    Author: Fernando Sumudu
    with small changes
    """
    import itertools
    inf_iters = tuple(itertools.cycle(enumerate(it)) for it in iters)
    num_iters = len(inf_iters)
    cur_val = [None]*num_iters

    first_v = True
    while True:
        i, p = 0, num_iters
        while p and not i:
            p -= 1
            i, cur_val[p] = next(inf_iters[p])

        if not p and not i:
            if first_v:
                first_v = False
            else:
                break

        yield cur_val


def sqrt_mod_iter(a, p, domain=int):
    """
    iterate over solutions to ``x**2 = a mod p``

    Parameters
    ==========

    a : integer
    p : positive integer
    domain : integer domain, ``int``, ``ZZ`` or ``Integer``

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import sqrt_mod_iter
    >>> list(sqrt_mod_iter(11, 43))
    [21, 22]
    """
    from sympy.polys.galoistools import gf_crt1, gf_crt2
    from sympy.polys.domains import ZZ
    a, p = as_int(a), abs(as_int(p))
    if isprime(p):
        a = a % p
        if a == 0:
            res = _sqrt_mod1(a, p, 1)
        else:
            res = _sqrt_mod_prime_power(a, p, 1)
        if res:
            if domain is ZZ:
                for x in res:
                    yield x
            else:
                for x in res:
                    yield domain(x)
    else:
        f = factorint(p)
        v = []
        pv = []
        for px, ex in f.items():
            if a % px == 0:
                rx = _sqrt_mod1(a, px, ex)
                if not rx:
                    raise StopIteration
            else:
                rx = _sqrt_mod_prime_power(a, px, ex)
                if not rx:
                    raise StopIteration
            v.append(rx)
            pv.append(px**ex)
        mm, e, s = gf_crt1(pv, ZZ)
        if domain is ZZ:
            for vx in _product(*v):
                r = gf_crt2(vx, pv, mm, e, s, ZZ)
                yield r
        else:
            for vx in _product(*v):
                r = gf_crt2(vx, pv, mm, e, s, ZZ)
                yield domain(r)


def _sqrt_mod_prime_power(a, p, k):
    """
    find the solutions to ``x**2 = a mod p**k`` when ``a % p != 0``

    Parameters
    ==========

    a : integer
    p : prime number
    k : positive integer

    References
    ==========

    [1] P. Hackman "Elementary Number Theory" (2009),  page 160
    [2] http://www.numbertheory.org/php/squareroot.html
    [3] [Gathen99]_

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _sqrt_mod_prime_power
    >>> _sqrt_mod_prime_power(11, 43, 1)
    [21, 22]
    """
    from sympy.core.numbers import igcdex
    from sympy.polys.domains import ZZ

    pk = p**k
    a = a % pk

    if k == 1:
        if p == 2:
            return [ZZ(a)]
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
            res = _sqrt_mod_tonelli_shanks(a, p)

        # ``_sqrt_mod_tonelli_shanks(a, p)`` is not deterministic;
        # sort to get always the same result
        return sorted([ZZ(res), ZZ(p - res)])

    if k > 1:
        f = factorint(a)
        # see Ref.[2]
        if p == 2:
            if a % 8 != 1:
                return None
            if k <= 3:
               s = set()
               for i in range(0, pk, 4):
                    s.add(1 + i)
                    s.add(-1 + i)
               return list(s)
            # according to Ref.[2] for k > 2 there are two solutions
            # (mod 2**k-1), that is four solutions (mod 2**k), which can be
            # obtained from the roots of x**2 = 0 (mod 8)
            rv = [ZZ(1), ZZ(3), ZZ(5), ZZ(7)]
            # hensel lift them to solutions of x**2 = 0 (mod 2**k)
            # if r**2 - a = 0 mod 2**nx but not mod 2**(nx+1)
            # then r + 2**(nx - 1) is a root mod 2**(nx+1)
            n = 3
            res = []
            for r in rv:
                nx = n
                while nx < k:
                    r1 = (r**2 - a) >> nx
                    if r1 % 2:
                        r = r + (1 << (nx - 1))
                    #assert (r**2 - a)% (1 << (nx + 1)) == 0
                    nx += 1
                if r not in res:
                    res.append(r)
                x = r + (1 << (k - 1))
                #assert (x**2 - a) % pk == 0
                if x < (1 << nx) and x not in res:
                    if (x**2 - a) % pk == 0:
                        res.append(x)
            return res
        rv = _sqrt_mod_prime_power(a, p, 1)
        if not rv:
            return None
        r = rv[0]
        fr = r**2 - a
        # hensel lifting with Newton iteration, see Ref.[3] chapter 9
        # with f(x) = x**2 - a; one has f'(a) != 0 (mod p) for p != 2
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
    find solution to ``x**2 == a mod p**n`` when ``a % p == 0``

    see http://www.numbertheory.org/php/squareroot.html
    """
    pn = p**n
    a = a % pn
    if a == 0:
        # case gcd(a, p**k) = p**n
        m = n // 2
        if n % 2 == 1:
            pm1 = p**(m + 1)
            def _iter0a():
                i = 0
                while i < pn:
                    yield i
                    i += pm1
            return _iter0a()
        else:
            pm = p**m
            def _iter0b():
                i = 0
                while i < pn:
                    yield i
                    i += pm
            return _iter0b()

    # case gcd(a, p**k) = p**r, r < n
    f = factorint(a)
    r = f[p]
    if r % 2 == 1:
        return None
    m = r // 2
    a1 = a >> r
    if p == 2:
        if n - r == 1:
            pnm1 = 1 << (n - m + 1)
            pm1 = 1 << (m + 1)
            def _iter1():
                k = 1 << (m + 2)
                i = 1 << m
                while i < pnm1:
                    j = i
                    while j < pn:
                        yield j
                        j += k
                    i += pm1
            return _iter1()
        if n - r == 2:
            res = _sqrt_mod_prime_power(a1, p, n - r)
            if res is None:
                return None
            pnm = 1 << (n - m)
            def _iter2():
                s = set()
                for r in res:
                    i = 0
                    while i < pn:
                        x = (r << m) + i
                        if x not in s:
                            s.add(x)
                            yield x
                        i += pnm
            return _iter2()
        if n - r > 2:
            res = _sqrt_mod_prime_power(a1, p, n - r)
            if res is None:
                return None
            pnm1 = 1 << (n - m - 1)
            def _iter3():
                s = set()
                for r in res:
                    i = 0
                    while i < pn:
                        x = ((r << m) + i) % pn
                        if x not in s:
                            s.add(x)
                            yield x
                        i += pnm1
            return _iter3()
    else:
        m = r // 2
        a1 = a // p**r
        res1 = _sqrt_mod_prime_power(a1, p, n - r)
        if res1 is None:
            return None
        pm = p**m
        pnr = p**(n-r)
        pnm = p**(n-m)

        def _iter4():
            s = set()
            pm = p**m
            for rx in res1:
                i = 0
                while i < pnm:
                    x = ((rx + i) % pn)
                    if x not in s:
                        s.add(x)
                        yield x*pm
                    i += pnr
        return _iter4()


def is_quad_residue(a, p):
    """
    Returns True if ``a`` (mod ``p``) is in the set of squares mod ``p``,
    i.e a % p in set([i**2 % p for i in range(p)]). If ``p`` is an odd
    prime, an iterative method is used to make the determination:

    >>> from sympy.ntheory import is_quad_residue
    >>> sorted(set([i**2 % 7 for i in range(7)]))
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
        r = sqrt_mod(a, p)
        if r is None:
            return False
        else:
            return True

    return pow(a, (p - 1) // 2, p) == 1


def is_nthpow_residue(a, n, m):
    """
    Returns True if ``x**n == a (mod m)`` has solutions.

    References
    ==========

    P. Hackman "Elementary Number Theory" (2009),  page 76
    """
    a, n, m = [as_int(i) for i in (a, n, m)]
    if m <= 0:
        raise ValueError('m must be > 0')
    if n < 0:
        raise ValueError('n must be >= 0')
    if a < 0:
        raise ValueError('a must be >= 0')
    if n == 0:
        if m == 1:
            return False
        return a == 1
    if n == 1:
        return True
    if n == 2:
        return is_quad_residue(a, m)
    return _is_nthpow_residue_bign(a, n, m)


def _is_nthpow_residue_bign(a, n, m):
    """Returns True if ``x**n == a (mod m)`` has solutions for n > 2."""
    # assert n > 2
    # assert a > 0 and m > 0
    if primitive_root(m) is None:
        # assert m >= 8
        for prime, power in factorint(m).items():
            if not _is_nthpow_residue_bign_prime_power(a, n, prime, power):
                return False
        return True
    f = totient(m)
    k = f // igcd(f, n)
    return pow(a, k, m) == 1


def _is_nthpow_residue_bign_prime_power(a, n, p, k):
    """Returns True/False if a solution for ``x**n == a (mod(p**k))``
    does/doesn't exist."""
    # assert a > 0
    # assert n > 2
    # assert p is prime
    # assert k > 0
    if a % p:
        if p != 2:
            return _is_nthpow_residue_bign(a, n, pow(p, k))
        if n & 1:
            return True
        c = trailing(n)
        return a % pow(2, min(c + 2, k)) == 1
    else:
        a %= pow(p, k)
        if not a:
            return True
        mu = multiplicity(p, a)
        if mu % n:
            return False
        pm = pow(p, mu)
        return _is_nthpow_residue_bign_prime_power(a//pm, n, p, k - mu)


def _nthroot_mod2(s, q, p):
    f = factorint(q)
    v = []
    for b, e in f.items():
        v.extend([b]*e)
    for qx in v:
        s = _nthroot_mod1(s, qx, p, False)
    return s


def _nthroot_mod1(s, q, p, all_roots):
    """
    Root of ``x**q = s mod p``, ``p`` prime and ``q`` divides ``p - 1``

    References
    ==========

    [1] A. M. Johnston "A Generalized qth Root Algorithm"
    """
    g = primitive_root(p)
    if not isprime(q):
        r = _nthroot_mod2(s, q, p)
    else:
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
        for t in range(p):
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
        return sqrt_mod(a, p , all_roots)
    f = totient(p)
    # see Hackman "Elementary Number Theory" (2009), page 76
    if not is_nthpow_residue(a, n, p):
        return None
    if not isprime(p) or primitive_root(p) == None:
        res = []
        moduli = []
        for prime, power in factorint(p).items():
            moduli.append(pow(prime, power))
            tmp = []
            _nthroot_mod_list(tmp, a, n, prime, power, all_roots)
            print (tmp)
            res.append(tmp)
        return sorted(_crt_cartesian(res, moduli))

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
        return sqrt_mod(a, p , all_roots)
    else:
        res = _nthroot_mod1(a, pa, p, all_roots)
    return res

def _nthroot_mod_list(roots, a, n, p, k, all_roots=False):
    """
    function return solution of x**n == a mod( p**k )

    """
    _roots = []
    if a % p:
        if p == 2:
            if k == 1:
                roots.append(1)
                return
            c = trailing(n)
            if k == 2:
                roots.append( a % 4 )
                if all_roots and c:
                    roots.append(3)
                return
            c = min(c, k - 2 )
            t = pow(2, k - 2)
            s = mod_inverse(r, t)
            if not c:
                roots.append(pow(a, s, pow(2, k) ))
                return
            t = a % pow(2, c + 2)
            root = 1
            pc = pow(2, c)
            pj = pc * 4
            for j in range(c+2, k):
                pj *= 2
                t = pow(root, pc, pj) - a
                if t % pj:
                    root += pow(2, j - c)
            pk = pow(2, k)
            root = pow(root, s, pk)
            if all_roots:
                t = pk / pc * root
                for i in range(2):
                    for j in range(pc):
                        roots.append(root)
                        root += t
                    root = t - root
            else:
                roots.append(root)
        else:
            #TODO implement _nthroot_mod_prime for general p**k
            #roots.extend(_nthroot_mod_prime(a, n, pow(p, k)) )
            return
    else:
        pk = pow(p, k)
        a %= pk
        m = 0
        pm = 0
        if not a:
            if not all_roots:
                roots.append(0)
                return
            _roots.append(0)
            if n >= k:
                m = k - 1
            else:
                m = k - 1 - ( k - 1 ) / n
            pm = pow(p, m)
        else:
            r = multiplicity(p, a)
            a /= pow(p, r)
            _nthroot_mod_list(_roots, a, n, p, k - r, all_roots)
            m = r / n
            pm = pow(p, m)
            if not all_roots:
                roots.append( _roots.back() * pow(p, m) )
                return
            _roots = [ pm*x for x in _roots]
            m = r - r / n
            pm = pow(p, m)

    pkm = pow(p, k - m )
    for x in _roots:
        root = x
        for i in range(pm):
            roots.append(root)
            root += pkm

def _crt_cartesian(rem, mod ):
    if len(mod) > len(rem):
        raiseerror
    if len(mod) == 0:
        raiseerror
    m = mod[0]
    R = rem[0]

    for i in range(1,len(mod)):
        rem2 = []
        s = mod_inverse(m, mod[i])
        _m = m
        m *= mod[i]
        for elem in R:
            for k in rem[i]:
                r = elem
                r += _m*s*(k - r)
                r %= m
                rem2.append(r)
        R = rem2
    return R

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
    for i in range(p // 2 + 1):
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
    >>> sorted(set([i**2 % 7 for i in range(7)]))
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


class mobius(Function):
    """
    Möbius function maps natural number to {-1, 0, 1}

    It is defined as follows:
        1) `1` if `n = 1`.
        2) `0` if `n` has a squared prime factor.
        3) `(-1)^k` if `n` is a square-free positive integer with `k`
           number of prime factors.

    It is an important multiplicative function in number theory
    and combinatorics.  It has applications in mathematical series,
    algebraic number theory and also physics (Fermion operator has very
    concrete realization with Möbius Function model).

    Parameters
    ==========

    n : positive integer

    Examples
    ========

    >>> from sympy.ntheory import mobius
    >>> mobius(13*7)
    1
    >>> mobius(1)
    1
    >>> mobius(13*7*5)
    -1
    >>> mobius(13**2)
    0

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/M%C3%B6bius_function
    .. [2] Thomas Koshy "Elementary Number Theory with Applications"

    """
    @classmethod
    def eval(cls, n):
        if n.is_integer:
            if n.is_positive is not True:
                raise ValueError("n should be a positive integer")
        else:
            raise TypeError("n should be an integer")
        if n.is_prime:
            return S.NegativeOne
        elif n is S.One:
            return S.One
        elif n.is_Integer:
            a = factorint(n)
            if any(i > 1 for i in a.values()):
                return S.Zero
            return S.NegativeOne**len(a)
