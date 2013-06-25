from sympy.ntheory import nextprime
from sympy.ntheory.modular import crt
from sympy.polys.galoistools import gf_gcd, gf_degree

def modgcd(f, g):
    """
    Computes the GCD of two polynomials in ``Z[x]`` using a modular
    algorithm.

    The algorithm computes the GCD of two univariate integer polynomials
    ``f`` and ``g`` by computing the GCD in ``Z_p[x]`` for suitable primes
    ``p`` and then reconstructing the coefficients with the Chinese
    Remainder Theorem. Trial division is only made for candidates which
    are very likely the desired GCD.

    Parameters
    ==========

    f : PolyElement
        univariate integer polynomial

    g : PolyElement
        univariate integer polynomial

    Returns
    =======

    h : PolyElement
        GCD of the polynomials f and g
    cff : PolyElement
        cofactor of f, i.e. quo(f, h)
    cfg : PolyElement
        cofactor of g, i.e. quo(g, h)

    Examples
    ========

    >>> from sympy.polys.modulargcd import modgcd
    >>> from sympy.polys import ring, ZZ

    >>> R, x = ring("x", ZZ)

    >>> f = x**5 - 1
    >>> g = x - 1

    >>> h, cff, cfg = modgcd(f, g)
    >>> h, cff, cfg
    (x - 1, x**4 + x**3 + x**2 + x + 1, 1)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    >>> f = 6*x**2 - 6
    >>> g = 2*x**2 + 4*x + 2

    >>> h, cff, cfg = modgcd(f, g)
    >>> h, cff, cfg
    (2*x + 2, 3*x - 3, x + 1)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    References
    ==========

    1. [Monagan00]_

    """
    assert f.ring == g.ring and f.ring.domain.is_ZZ

    result = _trivial_gcd(f, g)
    if result is not None:
        return result

    ring = f.ring

    cf, f = f.primitive()
    cg, g = g.primitive()
    ch = ring.domain.gcd(cf, cg)

    bound = _degree_bound(f, g)
    if bound == 0:
        return ch, f.mul_ground(cf/ch), g.mul_ground(cg/ch)

    gamma = ring.domain.gcd(f.LC, g.LC)
    m = 1
    p = 1

    while True:
        p = nextprime(p)
        while gamma % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        hp = _gf_gcd(fp, gp, p)
        (deghp,) = hp.LM # TODO: use hp.degree() instead

        if deghp > bound:
            continue
        elif deghp < bound:
            m = 1
            bound = deghp
            continue

        hp = hp * gamma
        hp = hp.trunc_ground(p)
        if m == 1:
            m = p
            hlastm = hp
            continue

        hm = _chinese_remainder_reconstruction(hp, hlastm, p, m)
        m *= p

        if not hm == hlastm:
            hlastm = hm
            continue

        h = hm.quo_ground(hm.content())
        if not f.rem(h) and not g.rem(h):
            if h.LC < 0:
                h = -h
            h = h.mul_ground(ch)
            cff = f.mul_ground(cf).quo(h)
            cfg = g.mul_ground(cg).quo(h)
            return h, cff, cfg


def _degree_bound(f, g):
    """
    Compute an upper bound for the degree of the GCD of two univariate
    integer polynomials ``f`` and ``g``.

    The function chooses a suitable prime ``p`` and computes the GCD of
    ``f`` and ``g`` in ``Z_p[x]``. The choice of ``p`` guarantees that
    the degree in ``Z_p[x]`` is greater than or equal to the degree in
    ``Z[x]``.

    Parameters
    ==========

    f : PolyElement
        univariate integer polynomial
    g : PolyElement
        univariate integer polynomial

    """
    gamma = f.ring.domain.gcd(f.LC, g.LC)
    p = 1

    while True:
        p = nextprime(p)
        while gamma % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        hp = _gf_gcd(fp, gp, p)
        (deghp,) = hp.LM # TODO: use hp.degree() instead
        return deghp


def _chinese_remainder_reconstruction(hp, hq, p, q):
    """
    Construct a polynomial ``hpq`` in ``Z_{p*q}[x]`` such that::

        hpq = hp (mod p)
        hpq = hq (mod q)

    for relatively prime integers ``p`` and ``q`` and polynomials
    ``hp`` and ``hq`` in ``Z_p[x]`` and ``Z_q[x]`` respectively.

    The coefficients of the polynomial ``hpq`` are computed with the
    Chinese Remainder Theorem. The symmetric representation in ``Z_p[x]``,
    ``Z_q[x]`` and ``Z_{p*q}[x]`` is used. It is assumed that ``hp`` and
    ``hq`` have the same degree.

    Parameters
    ==========

    hp : PolyElement
        univariate integer polynomial with coefficients in ``Z_p``
    hq : PolyElement
        univariate integer polynomial with coefficients in ``Z_q``
    p : Integer
        modulus of hp, relatively prime to q
    q : Integer
        modulus of hq, relatively prime to p

    Examples
    ========

    >>> from sympy.polys.modulargcd import _chinese_remainder_reconstruction
    >>> from sympy.polys import ring, ZZ

    >>> R, x = ring("x", ZZ)
    >>> p = 3
    >>> q = 5

    >>> hp = -x**3 - 1
    >>> hq = 2*x**3 - 2*x**2 + x

    >>> hpq = _chinese_remainder_reconstruction(hp, hq, p, q)
    >>> hpq
    2*x**3 + 3*x**2 + 6*x + 5

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    """
    (n,) = hp.LM # TODO: use hp.degree() instead
    x = hp.ring.gens[0]
    hpq = hp.ring.zero

    for i in range(n+1):
        hpq[(i,)] = crt([p, q], [hp.coeff(x**i), hq.coeff(x**i)], symmetric=True)[0]

    hpq.strip_zero()
    return hpq


# should I implement a GCD in GF(p) for PolyElement as well?
def _gf_gcd(fp, gp, p):
    """Compute the GCD of two univariate polynomials in ``Z_p[x]``."""
    ring = fp.ring
    densehp = gf_gcd(fp.to_dense(), gp.to_dense(), p, ring.domain)
    return ring.from_dense(densehp)


def _trivial_gcd(f, g):
    """
    Compute the GCD of two univariate polynomials for trivial cases,
    i.e. when one or both polynomials are zero.
    """
    ring = f.ring

    if not (f or g):
        return ring.zero, ring.zero, ring.zero
    elif not f:
        if g.LC < 0:
            return -g, ring.zero, -ring.one
        else:
            return g, ring.zero, ring.one
    elif not g:
        if f.LC < 0:
            return -f, -ring.one, ring.zero
        else:
            return f, ring.one, ring.zero
    return None
