from __future__ import annotations

from sympy.core.function import Function
from sympy.core.numbers import igcd, igcdex, mod_inverse
from sympy.core.power import isqrt
from sympy.core.singleton import S
from sympy.polys import Poly
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_crt1, gf_crt2, linear_congruence
from .primetest import isprime
from .factor_ import factorint, trailing, totient, multiplicity, perfect_power
from sympy.utilities.misc import as_int
from sympy.core.random import _randint, randint

from itertools import cycle, product


def n_order(a, n):
    """Returns the order of ``a`` modulo ``n``.

    The order of ``a`` modulo ``n`` is the smallest integer
    ``k`` such that ``a**k`` leaves a remainder of 1 with ``n``.

    Parameters
    ==========

    a : integer
    n : integer, n > 1. a and n should be relatively prime

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
    if n <= 1:
        raise ValueError("n should be an integer greater than 1")
    a = a % n
    # Trivial
    if a == 1:
        return 1
    if igcd(a, n) != 1:
        raise ValueError("The two numbers should be relatively prime")
    # We want to calculate
    # order = totient(n), factors = factorint(order)
    factors = defaultdict(int)
    for px, kx in factorint(n).items():
        if kx > 1:
            factors[px] += kx - 1
        for py, ky in factorint(px - 1).items():
            factors[py] += ky
    order = 1
    for px, kx in factors.items():
        order *= px**kx
    # Now the `order` is the order of the group.
    # The order of `a` divides the order of the group.
    for p, e in factors.items():
        for _ in range(e):
            if pow(a, order // p, n) == 1:
                order //= p
            else:
                break
    return order


def _primitive_root_prime_iter(p):
    """
    Generates the primitive roots for a prime ``p``

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _primitive_root_prime_iter
    >>> list(_primitive_root_prime_iter(19))
    [2, 3, 10, 13, 14, 15]

    References
    ==========

    .. [1] W. Stein "Elementary Number Theory" (2011), page 44

    """
    # it is assumed that p is an int
    v = [(p - 1) // i for i in factorint(p - 1).keys()]
    a = 2
    while a < p:
        for pw in v:
            # a TypeError below may indicate that p was not an int
            if pow(a, pw, p) == 1:
                break
        else:
            yield a
        a += 1


def primitive_root(p):
    """
    Returns the smallest primitive root or None.

    Parameters
    ==========

    p : positive integer

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import primitive_root
    >>> primitive_root(19)
    2

    References
    ==========

    .. [1] W. Stein "Elementary Number Theory" (2011), page 44
    .. [2] P. Hackman "Elementary Number Theory" (2009), Chapter C

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
    Returns True if ``a`` is a primitive root of ``p``.

    ``a`` is said to be the primitive root of ``p`` if gcd(a, p) == 1 and
    totient(p) is the smallest positive number s.t.

        a**totient(p) cong 1 mod(p)

    Parameters
    ==========

    a : integer
    p : integer, p > 1. a and p should be relatively prime

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
    if p <= 1:
        raise ValueError("p should be an integer greater than 1")
    a = a % p
    if igcd(a, p) != 1:
        raise ValueError("The two numbers should be relatively prime")
    # Primitive root of p exist only for
    # p = 2, 4, q**e, 2*q**e (q is odd prime)
    if p <= 4:
        # The primitive root is only p-1.
        return a == p - 1
    t = trailing(p)
    if t > 1:
        return False
    q = p >> t
    if isprime(q):
        group_order = q - 1
        factors = set(factorint(q - 1).keys())
    else:
        m = perfect_power(q)
        if not m:
            return False
        q, e = m
        if not isprime(q):
            return False
        group_order = q**(e - 1)*(q - 1)
        factors = set(factorint(q - 1).keys())
        factors.add(q)
    return all(pow(a, group_order // prime, p) != 1 for prime in factors)


def _sqrt_mod_tonelli_shanks(a, p):
    """
    Returns the square root in the case of ``p`` prime with ``p == 1 (mod 8)``

    References
    ==========

    .. [1] R. Crandall and C. Pomerance "Prime Numbers", 2nt Ed., page 101

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
    Find a root of ``x**2 = a mod p``.

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
        return sorted(sqrt_mod_iter(a, p))
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
    Cartesian product generator

    Notes
    =====

    Unlike itertools.product, it works also with iterables which do not fit
    in memory. See https://bugs.python.org/issue10109

    Author: Fernando Sumudu
    with small changes
    """
    inf_iters = tuple(cycle(enumerate(it)) for it in iters)
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
    Iterate over solutions to ``x**2 = a mod p``.

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
    a, p = as_int(a), abs(as_int(p))
    if isprime(p):
        a = a % p
        if a == 0:
            res = _sqrt_mod1(a, p, 1)
        else:
            res = _sqrt_mod_prime_power(a, p, 1)
        if res:
            if domain is ZZ:
                yield from res
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
                    return
            else:
                rx = _sqrt_mod_prime_power(a, px, ex)
                if not rx:
                    return
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
    Find the solutions to ``x**2 = a mod p**k`` when ``a % p != 0``

    Parameters
    ==========

    a : integer
    p : prime number
    k : positive integer

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _sqrt_mod_prime_power
    >>> _sqrt_mod_prime_power(11, 43, 1)
    [21, 22]

    References
    ==========

    .. [1] P. Hackman "Elementary Number Theory" (2009), page 160
    .. [2] http://www.numbertheory.org/php/squareroot.html
    .. [3] [Gathen99]_
    """
    pk = p**k
    a = a % pk

    if k == 1:
        if p == 2:
            return [ZZ(a)]
        if not (a % p < 2 or pow(a, (p - 1) // 2, p) == 1):
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
    Find solution to ``x**2 == a mod p**n`` when ``a % p == 0``

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
    i.e a % p in set([i**2 % p for i in range(p)]).

    Examples
    ========

    If ``p`` is an odd
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

    .. [1] P. Hackman "Elementary Number Theory" (2009), page 76

    """
    a = a % m
    a, n, m = as_int(a), as_int(n), as_int(m)
    if m <= 0:
        raise ValueError('m must be > 0')
    if n < 0:
        raise ValueError('n must be >= 0')
    if n == 0:
        if m == 1:
            return False
        return a == 1
    if a == 0:
        return True
    if n == 1:
        return True
    if n == 2:
        return is_quad_residue(a, m)
    return _is_nthpow_residue_bign(a, n, m)


def _is_nthpow_residue_bign(a, n, m):
    r"""Returns True if `x^n = a \pmod{n}` has solutions for `n > 2`."""
    # assert n > 2
    # assert a > 0 and m > 0
    if primitive_root(m) is None or igcd(a, m) != 1:
        # assert m >= 8
        for prime, power in factorint(m).items():
            if not _is_nthpow_residue_bign_prime_power(a, n, prime, power):
                return False
        return True
    f = totient(m)
    k = int(f // igcd(f, n))
    return pow(a, k, int(m)) == 1


def _is_nthpow_residue_bign_prime_power(a, n, p, k):
    r"""Returns True/False if a solution for `x^n = a \pmod{p^k}`
    does/does not exist."""
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

    .. [1] A. M. Johnston "A Generalized qth Root Algorithm"

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
        r1 = pow(s, x, p)
        s1 = pow(s, f, p)
        h = pow(g, f*q, p)
        t = discrete_log(p, s1, h)
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



def _help(m, prime_modulo_method, diff_method, expr_val):
    """
    Helper function for _nthroot_mod_composite and polynomial_congruence.

    Parameters
    ==========

    m : positive integer
    prime_modulo_method : function to calculate the root of the congruence
    equation for the prime divisors of m
    diff_method : function to calculate derivative of expression at any
    given point
    expr_val : function to calculate value of the expression at any
    given point
    """
    from sympy.ntheory.modular import crt
    f = factorint(m)
    dd = {}
    for p, e in f.items():
        tot_roots = set()
        if e == 1:
            tot_roots.update(prime_modulo_method(p))
        else:
            for root in prime_modulo_method(p):
                diff = diff_method(root, p)
                if diff != 0:
                    ppow = p
                    m_inv = mod_inverse(diff, p)
                    for j in range(1, e):
                        ppow *= p
                        root = (root - expr_val(root, ppow) * m_inv) % ppow
                    tot_roots.add(root)
                else:
                    new_base = p
                    roots_in_base = {root}
                    while new_base < pow(p, e):
                        new_base *= p
                        new_roots = set()
                        for k in roots_in_base:
                            if expr_val(k, new_base)!= 0:
                                continue
                            while k not in new_roots:
                                new_roots.add(k)
                                k = (k + (new_base // p)) % new_base
                        roots_in_base = new_roots
                    tot_roots = tot_roots | roots_in_base
        if tot_roots == set():
            return []
        dd[pow(p, e)] = tot_roots
    a = []
    m = []
    for x, y in dd.items():
        m.append(x)
        a.append(list(y))
    return sorted({crt(m, list(i))[0] for i in product(*a)})


def _nthroot_mod_composite(a, n, m):
    """
    Find the solutions to ``x**n = a mod m`` when m is not prime.
    """
    return _help(m,
        lambda p: nthroot_mod(a, n, p, True),
        lambda root, p: (pow(root, n - 1, p) * (n % p)) % p,
        lambda root, p: (pow(root, n, p) - a) % p)


def nthroot_mod(a, n, p, all_roots=False):
    """
    Find the solutions to ``x**n = a mod p``.

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
    a = a % p
    a, n, p = as_int(a), as_int(n), as_int(p)

    if n == 2:
        return sqrt_mod(a, p, all_roots)
    # see Hackman "Elementary Number Theory" (2009), page 76
    if not isprime(p):
        return _nthroot_mod_composite(a, n, p)
    if a % p == 0:
        return [0]
    if not is_nthpow_residue(a, n, p):
        return [] if all_roots else None
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
        return sqrt_mod(a, p, all_roots)
    else:
        res = _nthroot_mod1(a, pa, p, all_roots)
    return res


def quadratic_residues(p) -> list[int]:
    """
    Returns the list of quadratic residues.

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import quadratic_residues
    >>> quadratic_residues(7)
    [0, 1, 2, 4]
    """
    p = as_int(p)
    r = {pow(i, 2, p) for i in range(p // 2 + 1)}
    return sorted(r)


def legendre_symbol(a, p):
    r"""
    Returns the Legendre symbol `(a / p)`.

    For an integer ``a`` and an odd prime ``p``, the Legendre symbol is
    defined as

    .. math ::
        \genfrac(){}{}{a}{p} = \begin{cases}
             0 & \text{if } p \text{ divides } a\\
             1 & \text{if } a \text{ is a quadratic residue modulo } p\\
            -1 & \text{if } a \text{ is a quadratic nonresidue modulo } p
        \end{cases}

    Parameters
    ==========

    a : integer
    p : odd prime

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
    if pow(a, (p - 1) // 2, p) == 1:
        return 1
    return -1


def jacobi_symbol(m, n):
    r"""
    Returns the Jacobi symbol `(m / n)`.

    For any integer ``m`` and any positive odd integer ``n`` the Jacobi symbol
    is defined as the product of the Legendre symbols corresponding to the
    prime factors of ``n``:

    .. math ::
        \genfrac(){}{}{m}{n} =
            \genfrac(){}{}{m}{p^{1}}^{\alpha_1}
            \genfrac(){}{}{m}{p^{2}}^{\alpha_2}
            ...
            \genfrac(){}{}{m}{p^{k}}^{\alpha_k}
            \text{ where } n =
                p_1^{\alpha_1}
                p_2^{\alpha_2}
                ...
                p_k^{\alpha_k}

    Like the Legendre symbol, if the Jacobi symbol `\genfrac(){}{}{m}{n} = -1`
    then ``m`` is a quadratic nonresidue modulo ``n``.

    But, unlike the Legendre symbol, if the Jacobi symbol
    `\genfrac(){}{}{m}{n} = 1` then ``m`` may or may not be a quadratic residue
    modulo ``n``.

    Parameters
    ==========

    m : integer
    n : odd positive integer

    Examples
    ========

    >>> from sympy.ntheory import jacobi_symbol, legendre_symbol
    >>> from sympy import S
    >>> jacobi_symbol(45, 77)
    -1
    >>> jacobi_symbol(60, 121)
    1

    The relationship between the ``jacobi_symbol`` and ``legendre_symbol`` can
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
    if n < 0 or not n % 2:
        raise ValueError("n should be an odd positive integer")
    if m < 0 or m > n:
        m %= n
    if not m:
        return int(n == 1)
    if n == 1 or m == 1:
        return 1
    if igcd(m, n) != 1:
        return 0

    j = 1
    while m != 0:
        while m % 2 == 0 and m > 0:
            m >>= 1
            if n % 8 in [3, 5]:
                j = -j
        m, n = n, m
        if m % 4 == n % 4 == 3:
            j = -j
        m %= n
    return j


class mobius(Function):
    """
    Mobius function maps natural number to {-1, 0, 1}

    It is defined as follows:
        1) `1` if `n = 1`.
        2) `0` if `n` has a squared prime factor.
        3) `(-1)^k` if `n` is a square-free positive integer with `k`
           number of prime factors.

    It is an important multiplicative function in number theory
    and combinatorics.  It has applications in mathematical series,
    algebraic number theory and also physics (Fermion operator has very
    concrete realization with Mobius Function model).

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

    .. [1] https://en.wikipedia.org/wiki/M%C3%B6bius_function
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


def _discrete_log_trial_mul(n, a, b, order=None):
    """
    Trial multiplication algorithm for computing the discrete logarithm of
    ``a`` to the base ``b`` modulo ``n``.

    The algorithm finds the discrete logarithm using exhaustive search. This
    naive method is used as fallback algorithm of ``discrete_log`` when the
    group order is very small.

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _discrete_log_trial_mul
    >>> _discrete_log_trial_mul(41, 15, 7)
    3

    See Also
    ========

    discrete_log

    References
    ==========

    .. [1] "Handbook of applied cryptography", Menezes, A. J., Van, O. P. C., &
        Vanstone, S. A. (1997).
    """
    a %= n
    b %= n
    if order is None:
        order = n
    x = 1
    for i in range(order):
        if x == a:
            return i
        x = x * b % n
    raise ValueError("Log does not exist")


def _discrete_log_shanks_steps(n, a, b, order=None):
    """
    Baby-step giant-step algorithm for computing the discrete logarithm of
    ``a`` to the base ``b`` modulo ``n``.

    The algorithm is a time-memory trade-off of the method of exhaustive
    search. It uses `O(sqrt(m))` memory, where `m` is the group order.

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _discrete_log_shanks_steps
    >>> _discrete_log_shanks_steps(41, 15, 7)
    3

    See Also
    ========

    discrete_log

    References
    ==========

    .. [1] "Handbook of applied cryptography", Menezes, A. J., Van, O. P. C., &
        Vanstone, S. A. (1997).
    """
    a %= n
    b %= n
    if order is None:
        order = n_order(b, n)
    m = isqrt(order) + 1
    T = {}
    x = 1
    for i in range(m):
        T[x] = i
        x = x * b % n
    z = mod_inverse(b, n)
    z = pow(z, m, n)
    x = a
    for i in range(m):
        if x in T:
            return i * m + T[x]
        x = x * z % n
    raise ValueError("Log does not exist")


def _discrete_log_pollard_rho(n, a, b, order=None, retries=10, rseed=None):
    """
    Pollard's Rho algorithm for computing the discrete logarithm of ``a`` to
    the base ``b`` modulo ``n``.

    It is a randomized algorithm with the same expected running time as
    ``_discrete_log_shanks_steps``, but requires a negligible amount of memory.

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _discrete_log_pollard_rho
    >>> _discrete_log_pollard_rho(227, 3**7, 3)
    7

    See Also
    ========

    discrete_log

    References
    ==========

    .. [1] "Handbook of applied cryptography", Menezes, A. J., Van, O. P. C., &
        Vanstone, S. A. (1997).
    """
    a %= n
    b %= n

    if order is None:
        order = n_order(b, n)
    randint = _randint(rseed)

    for i in range(retries):
        aa = randint(1, order - 1)
        ba = randint(1, order - 1)
        xa = pow(b, aa, n) * pow(a, ba, n) % n

        c = xa % 3
        if c == 0:
            xb = a * xa % n
            ab = aa
            bb = (ba + 1) % order
        elif c == 1:
            xb = xa * xa % n
            ab = (aa + aa) % order
            bb = (ba + ba) % order
        else:
            xb = b * xa % n
            ab = (aa + 1) % order
            bb = ba

        for j in range(order):
            c = xa % 3
            if c == 0:
                xa = a * xa % n
                ba = (ba + 1) % order
            elif c == 1:
                xa = xa * xa % n
                aa = (aa + aa) % order
                ba = (ba + ba) % order
            else:
                xa = b * xa % n
                aa = (aa + 1) % order

            c = xb % 3
            if c == 0:
                xb = a * xb % n
                bb = (bb + 1) % order
            elif c == 1:
                xb = xb * xb % n
                ab = (ab + ab) % order
                bb = (bb + bb) % order
            else:
                xb = b * xb % n
                ab = (ab + 1) % order

            c = xb % 3
            if c == 0:
                xb = a * xb % n
                bb = (bb + 1) % order
            elif c == 1:
                xb = xb * xb % n
                ab = (ab + ab) % order
                bb = (bb + bb) % order
            else:
                xb = b * xb % n
                ab = (ab + 1) % order

            if xa == xb:
                r = (ba - bb) % order
                try:
                    e = mod_inverse(r, order) * (ab - aa) % order
                    if (pow(b, e, n) - a) % n == 0:
                        return e
                except ValueError:
                    pass
                break
    raise ValueError("Pollard's Rho failed to find logarithm")


def _discrete_log_pohlig_hellman(n, a, b, order=None):
    """
    Pohlig-Hellman algorithm for computing the discrete logarithm of ``a`` to
    the base ``b`` modulo ``n``.

    In order to compute the discrete logarithm, the algorithm takes advantage
    of the factorization of the group order. It is more efficient when the
    group order factors into many small primes.

    Examples
    ========

    >>> from sympy.ntheory.residue_ntheory import _discrete_log_pohlig_hellman
    >>> _discrete_log_pohlig_hellman(251, 210, 71)
    197

    See Also
    ========

    discrete_log

    References
    ==========

    .. [1] "Handbook of applied cryptography", Menezes, A. J., Van, O. P. C., &
        Vanstone, S. A. (1997).
    """
    from .modular import crt
    a %= n
    b %= n

    if order is None:
        order = n_order(b, n)

    f = factorint(order)
    l = [0] * len(f)

    for i, (pi, ri) in enumerate(f.items()):
        for j in range(ri):
            gj = pow(b, l[i], n)
            aj = pow(a * mod_inverse(gj, n), order // pi**(j + 1), n)
            bj = pow(b, order // pi, n)
            cj = discrete_log(n, aj, bj, pi, True)
            l[i] += cj * pi**j

    d, _ = crt([pi**ri for pi, ri in f.items()], l)
    return d


def discrete_log(n, a, b, order=None, prime_order=None):
    """
    Compute the discrete logarithm of ``a`` to the base ``b`` modulo ``n``.

    This is a recursive function to reduce the discrete logarithm problem in
    cyclic groups of composite order to the problem in cyclic groups of prime
    order.

    It employs different algorithms depending on the problem (subgroup order
    size, prime order or not):

        * Trial multiplication
        * Baby-step giant-step
        * Pollard's Rho
        * Pohlig-Hellman

    Examples
    ========

    >>> from sympy.ntheory import discrete_log
    >>> discrete_log(41, 15, 7)
    3

    References
    ==========

    .. [1] https://mathworld.wolfram.com/DiscreteLogarithm.html
    .. [2] "Handbook of applied cryptography", Menezes, A. J., Van, O. P. C., &
        Vanstone, S. A. (1997).

    """
    n, a, b = as_int(n), as_int(a), as_int(b)
    if order is None:
        order = n_order(b, n)

    if prime_order is None:
        prime_order = isprime(order)

    if order < 1000:
        return _discrete_log_trial_mul(n, a, b, order)
    elif prime_order:
        if order < 1000000000000:
            return _discrete_log_shanks_steps(n, a, b, order)
        return _discrete_log_pollard_rho(n, a, b, order)

    return _discrete_log_pohlig_hellman(n, a, b, order)



def quadratic_congruence(a, b, c, p):
    """
    Find the solutions to ``a x**2 + b x + c = 0 mod p.

    Parameters
    ==========

    a : int
    b : int
    c : int
    p : int
        A positive integer.
    """
    a = as_int(a)
    b = as_int(b)
    c = as_int(c)
    p = as_int(p)
    a = a % p
    b = b % p
    c = c % p

    if a == 0:
        return linear_congruence(b, -c, p)
    if p == 2:
        roots = []
        if c % 2 == 0:
            roots.append(0)
        if (a + b + c) % 2 == 0:
            roots.append(1)
        return roots
    if isprime(p):
        inv_a = mod_inverse(a, p)
        b *= inv_a
        c *= inv_a
        if b % 2 == 1:
            b = b + p
        d = ((b * b) // 4 - c) % p
        y = sqrt_mod(d, p, all_roots=True)
        res = set()
        for i in y:
            res.add((i - b // 2) % p)
        return sorted(res)
    y = sqrt_mod(b * b - 4 * a * c, 4 * a * p, all_roots=True)
    res = set()
    for i in y:
        root = linear_congruence(2 * a, i - b, 4 * a * p)
        for j in root:
            res.add(j % p)
    return sorted(res)


def _polynomial_congruence_prime(coefficients, p):
    """A helper function used by polynomial_congruence.
    It returns the root of a polynomial modulo prime number
    by naive search from [0, p).

    Parameters
    ==========

    coefficients : list of integers
    p : prime number
    """

    roots = []
    rank = len(coefficients)
    for i in range(0, p):
        f_val = 0
        for coeff in range(0,rank - 1):
            f_val = (f_val + pow(i, int(rank - coeff - 1), p) * coefficients[coeff]) % p
        f_val = f_val + coefficients[-1]
        if f_val % p == 0:
            roots.append(i)
    return roots


def _diff_poly(root, coefficients, p):
    """A helper function used by polynomial_congruence.
    It returns the derivative of the polynomial evaluated at the
    root (mod p).

    Parameters
    ==========

    coefficients : list of integers
    p : prime number
    root : integer
    """

    diff = 0
    rank = len(coefficients)
    for coeff in range(0, rank - 1):
        if not coefficients[coeff]:
            continue
        diff = (diff + pow(root, rank - coeff - 2, p)*(rank - coeff - 1)*
            coefficients[coeff]) % p
    return diff % p


def _val_poly(root, coefficients, p):
    """A helper function used by polynomial_congruence.
    It returns value of the polynomial at root (mod p).

    Parameters
    ==========

    coefficients : list of integers
    p : prime number
    root : integer
    """
    rank = len(coefficients)
    f_val = 0
    for coeff in range(0, rank - 1):
        f_val = (f_val + pow(root, rank - coeff - 1, p)*
            coefficients[coeff]) % p
    f_val = f_val + coefficients[-1]
    return f_val % p


def _valid_expr(expr):
    """
    return coefficients of expr if it is a univariate polynomial
    with integer coefficients else raise a ValueError.
    """

    if not expr.is_polynomial():
        raise ValueError("The expression should be a polynomial")
    polynomial = Poly(expr)
    if not polynomial.is_univariate:
        raise ValueError("The expression should be univariate")
    if not polynomial.domain == ZZ:
        raise ValueError("The expression should should have integer coefficients")
    return polynomial.all_coeffs()


def polynomial_congruence(expr, m):
    """
    Find the solutions to a polynomial congruence equation modulo m.

    Parameters
    ==========

    coefficients : Coefficients of the Polynomial
    m : positive integer

    Examples
    ========

    >>> from sympy.ntheory import polynomial_congruence
    >>> from sympy.abc import x
    >>> expr = x**6 - 2*x**5 -35
    >>> polynomial_congruence(expr, 6125)
    [3257]
    """
    coefficients = _valid_expr(expr)
    coefficients = [num % m for num in coefficients]
    rank = len(coefficients)
    if rank == 3:
        return quadratic_congruence(*coefficients, m)
    if rank == 2:
        return quadratic_congruence(0, *coefficients, m)
    if coefficients[0] == 1 and 1 + coefficients[-1] == sum(coefficients):
        return nthroot_mod(-coefficients[-1], rank - 1, m, True)
    if isprime(m):
        return _polynomial_congruence_prime(coefficients, m)
    return _help(m,
        lambda p: _polynomial_congruence_prime(coefficients, p),
        lambda root, p: _diff_poly(root, coefficients, p),
        lambda root, p: _val_poly(root, coefficients, p))
