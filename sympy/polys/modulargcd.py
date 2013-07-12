from sympy.ntheory import nextprime
from sympy.ntheory.modular import crt
from sympy.polys.galoistools import gf_gcd


def _trivial_gcd(f, g):
    """
    Compute the GCD of two polynomials in trivial cases, i.e. when one
    or both polynomials are zero.
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


# TODO: avoid conversion to dense representation
def _gf_gcd(fp, gp, p):
    """
    Compute the GCD of two univariate polynomials in `\mathbb{Z}_p[x]`.
    """
    ring = fp.ring
    densehp = gf_gcd(fp.to_dense(), gp.to_dense(), p, ring.domain)
    return ring.from_dense(densehp)


def _degree_bound_univariate(f, g):
    """
    Compute an upper bound for the degree of the GCD of two univariate
    integer polynomials `f` and `g`.

    The function chooses a suitable prime `p` and computes the GCD of
    `f` and `g` in `\mathbb{Z}_p[x]`. The choice of `p` guarantees that
    the degree in `\mathbb{Z}_p[x]` is greater than or equal to the degree
    in `\mathbb{Z}[x]`.

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
        deghp = hp.degree()
        return deghp


def _chinese_remainder_reconstruction_univariate(hp, hq, p, q):
    """
    Construct a polynomial `h_{pq}` in `\mathbb{Z}_{p q}[x]` such that

    .. math ::

        h_{pq} = h_p \; \mathrm{mod} \, p

        h_{pq} = h_q \; \mathrm{mod} \, q

    for relatively prime integers `p` and `q` and polynomials
    `h_p` and `h_q` in `\mathbb{Z}_p[x]` and `\mathbb{Z}_q[x]`
    respectively.

    The coefficients of the polynomial `h_{pq}` are computed with the
    Chinese Remainder Theorem. The symmetric representation in
    `\mathbb{Z}_p[x]`, `\mathbb{Z}_q[x]` and `\mathbb{Z}_{p q}[x]` is used.
    It is assumed that `h_p` and `h_q` have the same degree.

    Parameters
    ==========

    hp : PolyElement
        univariate integer polynomial with coefficients in `\mathbb{Z}_p`
    hq : PolyElement
        univariate integer polynomial with coefficients in `\mathbb{Z}_q`
    p : Integer
        modulus of `h_p`, relatively prime to `q`
    q : Integer
        modulus of `h_q`, relatively prime to `p`

    Examples
    ========

    >>> from sympy.polys.modulargcd import _chinese_remainder_reconstruction_univariate
    >>> from sympy.polys import ring, ZZ

    >>> R, x = ring("x", ZZ)
    >>> p = 3
    >>> q = 5

    >>> hp = -x**3 - 1
    >>> hq = 2*x**3 - 2*x**2 + x

    >>> hpq = _chinese_remainder_reconstruction_univariate(hp, hq, p, q)
    >>> hpq
    2*x**3 + 3*x**2 + 6*x + 5

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    """
    n = hp.degree()
    x = hp.ring.gens[0]
    hpq = hp.ring.zero

    for i in xrange(n+1):
        hpq[(i,)] = crt([p, q], [hp.coeff(x**i), hq.coeff(x**i)], symmetric=True)[0]

    hpq.strip_zero()
    return hpq


def modgcd_univariate(f, g):
    """
    Computes the GCD of two polynomials in `\mathbb{Z}[x]` using a modular
    algorithm.

    The algorithm computes the GCD of two univariate integer polynomials
    `f` and `g` by computing the GCD in `\mathbb{Z}_p[x]` for suitable
    primes `p` and then reconstructing the coefficients with the Chinese
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
        GCD of the polynomials `f` and `g`
    cff : PolyElement
        cofactor of `f`, i.e. `\\frac{f}{h}`
    cfg : PolyElement
        cofactor of `g`, i.e. `\\frac{g}{h}`

    Examples
    ========

    >>> from sympy.polys.modulargcd import modgcd_univariate
    >>> from sympy.polys import ring, ZZ

    >>> R, x = ring("x", ZZ)

    >>> f = x**5 - 1
    >>> g = x - 1

    >>> h, cff, cfg = modgcd_univariate(f, g)
    >>> h, cff, cfg
    (x - 1, x**4 + x**3 + x**2 + x + 1, 1)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    >>> f = 6*x**2 - 6
    >>> g = 2*x**2 + 4*x + 2

    >>> h, cff, cfg = modgcd_univariate(f, g)
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

    bound = _degree_bound_univariate(f, g)
    if bound == 0:
        return ring(ch), f.mul_ground(cf // ch), g.mul_ground(cg // ch)

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
        deghp = hp.degree()

        if deghp > bound:
            continue
        elif deghp < bound:
            m = 1
            bound = deghp
            continue

        hp = hp.mul_ground(gamma).trunc_ground(p)
        if m == 1:
            m = p
            hlastm = hp
            continue

        hm = _chinese_remainder_reconstruction_univariate(hp, hlastm, p, m)
        m *= p

        if not hm == hlastm:
            hlastm = hm
            continue

        h = hm.quo_ground(hm.content())
        fquo, frem = f.div(h)
        gquo, grem = g.div(h)
        if not frem and not grem:
            if h.LC < 0:
                ch = -ch
            h = h.mul_ground(ch)
            cff = fquo.mul_ground(cf // ch)
            cfg = gquo.mul_ground(cg // ch)
            return h, cff, cfg
