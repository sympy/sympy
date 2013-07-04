from sympy.ntheory import nextprime
from sympy.ntheory.modular import crt
from sympy.polys.galoistools import gf_gcd, gf_from_dict
from sympy.polys.densebasic import dmp_swap
from sympy.polys.polytools import invert
from sympy.polys.polyerrors import ModularGCDFailed

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
        (deghp,) = hp.LM # TODO: use hp.degree() instead

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

        hm = _chinese_remainder_reconstruction(hp, hlastm, p, m)
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


def modgcd_bivariate(f, g):
    assert f.ring == g.ring and f.ring.domain.is_ZZ

    result = _trivial_gcd(f, g)
    if result is not None:
        return result

    ring = f.ring

    cf, f = f.primitive()
    cg, g = g.primitive()
    ch = ring.domain.gcd(cf, cg)

    xbound, ycontbound = _degree_bound_bivariate(f, g)
    if xbound == ycontbound == 0:
        return ring(ch), f.mul_ground(cf // ch), g.mul_ground(cg // ch)

    fswap = _swap(f, 1)
    gswap = _swap(g, 1)
    (degyf, _) = fswap.LM
    (degyg, _) = gswap.LM

    ybound, xcontbound = _degree_bound_bivariate(fswap, gswap)
    if ybound == xcontbound == 0:
        return ring(ch), f.mul_ground(cf // ch), g.mul_ground(cg // ch)

    #TODO: CHOOSE MAIN VARIABLE x HERE

    gamma1 = ring.domain.gcd(f.LC, g.LC)
    gamma2 = ring.domain.gcd(fswap.LC, gswap.LC)
    badprimes = gamma1 * gamma2
    m = 1
    p = 1

    while True:
        p = nextprime(p)
        while badprimes % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        contfp, fp = ring.dmp_primitive(fp)
        contgp, gp = ring.dmp_primitive(gp)
        conthp = _gf_gcd(contfp, contgp, p) # monic polynomial in Z_p[y]
        (degconthp,) = conthp.LM # TODO: use conthp.degree() instead

        if degconthp > ycontbound:
            continue
        elif degconthp < ycontbound:
            m = 1
            ycontbound = degconthp
            continue

        delta = _gf_gcd(ring.dmp_LC(fp), ring.dmp_LC(gp), p) # polynomial in Z_p[y]
        y = delta.ring.gens[0]
        (degcontfp,) = contfp.LM
        (degcontgp,) = contgp.LM
        (degdelta,) = delta.LM

        N = min(degyf - degcontfp, degyg - degcontgp,
            ybound - ycontbound + degdelta) + 1

        if p < N:
            continue

        n = 0
        evalpoints = []
        hpeval = []
        unlucky = False

        for a in range(p):
            deltaa = delta.evaluate(y, a)
            if not deltaa % p:
                continue

            y_ = ring.gens[1] # problem: y != y_
            fpa = fp.evaluate(y_, a).trunc_ground(p)
            gpa = gp.evaluate(y_, a).trunc_ground(p)
            hpa = _gf_gcd(fpa, gpa, p) # monic polynomial in Z_p[x]
            (deghpa,) = hpa.LM # TODO: use hpa.degree() instead

            if deghpa > xbound:
                continue
            elif deghpa < xbound:
                m = 1
                xbound = deghpa
                unlucky = True
                break

            hpa = hpa.mul_ground(deltaa).trunc_ground(p)
            evalpoints.append(a)
            hpeval.append(hpa)
            n += 1

            if n == N:
                break

        if unlucky:
            continue
        if n < N:
            continue

        hp = _interpolate_multivariate(evalpoints, hpeval, ring, p)

        hp = ring.dmp_primitive(hp)[1]
        hp = hp * conthp.set_ring(ring)
        degyhp = ring.dmp_degree_in(hp, 1)

        if degyhp > ybound:
            continue
        if degyhp < ybound:
            m = 1
            ybound = degyhp
            continue

        hp = hp.mul_ground(gamma1).trunc_ground(p)
        if m == 1:
            m = p
            hlastm = hp
            continue

        hm = _chinese_remainder_reconstruction_multivariate(hp, hlastm, p, m)
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


def _swap(f, i):
    ring = f.ring
    k = ring.ngens
    f = ring.from_dense(dmp_swap(f.to_dense(), 0, i, k-1, ring.domain))
    return f


def _degree_bound_bivariate(f, g):
    ring = f.ring

    gamma1 = ring.domain.gcd(f.LC, g.LC)
    gamma2 = ring.domain.gcd(fswap.LC, gswap.LC)
    badprimes = gamma1 * gamma2
    p = 1

    while True:
        p = nextprime(p)
        while badprimes % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        contfp, fp = ring.dmp_primitive(fp)
        contgp, gp = ring.dmp_primitive(gp)
        conthp = _gf_gcd(contfp, contgp, p) # polynomial in Z[y]
        (ycontbound,) = conthp.LM # TODO: use conthp.degree() instead

        delta = _gf_gcd(ring.dmp_LC(fp), ring.dmp_LC(gp), p) # polynomial in Z[y]
        y = delta.ring.gens[0]

        for a in range(p):
            if not delta.evaluate(y, a) % p:
                continue
            y_ = ring.gens[1] # problem: y != y_
            fpa = fp.evaluate(y_, a).trunc_ground(p)
            gpa = gp.evaluate(y_, a).trunc_ground(p)
            hpa = _gf_gcd(fpa, gpa, p)
            (xbound,) = hpa.LM # TODO: use hpa.degree() instead
            return xbound, ycontbound

        return min(fp.LM[0], gp.LM[0]), ycontbound


def _chinese_remainder_reconstruction_multivariate(hp, hq, p, q):
    hpmonoms = set(hp.monoms())
    hqmonoms = set(hq.monoms())
    monoms = hpmonoms.intersection(hqmonoms)
    hpmonoms.difference_update(monoms)
    hqmonoms.difference_update(monoms)

    zero = hp.ring.domain.zero

    hpq = hp.ring.zero

    for monom in monoms:
        hpq[monom] = crt([p, q], [hp[monom], hq[monom]], symmetric=True)[0]
    for monom in hpmonoms:
        hpq[monom] = crt([p, q], [hp[monom], zero], symmetric=True)[0]
    for monom in hqmonoms:
        hpq[monom] = crt([p, q], [zero, hq[monom]], symmetric=True)[0]

    return hpq


def _interpolate_multivariate(evalpoints, hpeval, ring, p):
    hp = ring.zero
    k = ring.ngens
    y = ring.gens[k-1]
    for a, hpa in zip(evalpoints, hpeval):
        numer = ring.one
        denom = ring.domain.one
        for b in evalpoints:
            if b == a:
                continue

            numer *= y - b
            denom *= a - b

        coeff = numer.mul_ground(invert(denom, p))
        hp += hpa.set_ring(ring) * coeff

    return hp.trunc_ground(p)


# f, g polynomials in Z_p[x0, ..., x{k-2}, y] ~ Z_p[y][x0, ..., x{k-2}]
# X = (x0, ..., x{k-2})

def brown_p(f, g, p):
    """
    Compute the GCD and cofactors of two polynomials in
    Z_p[x0, ..., x{k-1}].
    """
    ring = f.ring
    k = ring.ngens

    if k == 1:
        h = _gf_gcd(f, g, p).trunc_ground(p)
        cff = f.quo(h)
        cfg = g.quo(h)
        return h, cff, cfg

    contf, f = _primitive(f, p)
    contg, g = _primitive(g, p)
    conth = _gf_gcd(contf, contg, p) # polynomial in Z_p[y]
    degyf = ring.dmp_degree_in(f, k-1)
    degyg = ring.dmp_degree_in(g, k-1)

    lcf = _LC(f)
    lcg = _LC(g)
    gamma = _gf_gcd(lcf, lcg, p) # polynomial in Z_p[y]
    (deggamma,) = gamma.LM
    y = gamma.ring.gens[0]

    Xbound = min(_deg(f), _deg(g))
    zero = tuple(0 for i in range(k-1))
    nuf = deggamma + degyf
    nug = deggamma + degyg
    N = max(nuf, nug) + 1
    n = 0
    evalpoints = []
    heval = []
    cffeval = []
    cfgeval = []

    for a in range(p):
        if not lcf.evaluate(y, a) * lcg.evaluate(y, a) % p:
            continue

        gammaa = gamma.evaluate(y, a) % p
        y1 = ring.gens[k-1] # problem: y != y1
        fa = f.evaluate(y1, a).mul_ground(gammaa).trunc_ground(p)
        ga = g.evaluate(y1, a).mul_ground(gammaa).trunc_ground(p)

        res = brown_p(fa, ga, p) # polynomials in Z_p[X], symmetric range
        if res == 'algorithm failed':
            continue

        ha, cffa, cfga = res
        degh = ha.LM

        if degh == zero:
            h = conth.set_ring(ring)
            # problem: cannot multiply a polynomial in Z_p[y] to a polynomial in Z_p[X,y],
            # i guess it is because the generators 'y' are not the same
            cff = contf.quo(conth).set_ring(ring) * f
            cfg = contg.quo(conth).set_ring(ring) * g
            return h.trunc_ground(p), cff.trunc_ground(p), cfg.trunc_ground(p)

        elif degh > Xbound:
            continue
        elif degh < Xbound:
            n = 0
            Xbound = degh

        evalpoints.append(a)
        heval.append(ha)
        cffeval.append(cffa)
        cfgeval.append(cfga)
        n += 1

        if n == N:
            h = _interpolate_multivariate(evalpoints, heval, ring, p)
            cff = _interpolate_multivariate(evalpoints, cffeval, ring, p)
            cfg = _interpolate_multivariate(evalpoints, cfgeval, ring, p)
            degyh = ring.dmp_degree_in(h, k-1)
            degycff = ring.dmp_degree_in(cff, k-1)
            degycfg = ring.dmp_degree_in(cfg, k-1)

            if degyh + degycff != nuf or degyh + degycfg != nug:
                n = 0
                continue

            h = _primitive(h, p)[1]
            lch = _LC(h).set_ring(ring)
            cff = contf.quo(conth).set_ring(ring) * cff.quo(lch)
            cfg = contg.quo(conth).set_ring(ring) * cfg.quo(lch)

            return conth.set_ring(ring) * h, cff.trunc_ground(p), cfg.trunc_ground(p)

    return 'algorithm failed'


def _primitive(f, p):
    """
    Compute the content and the primitive part of a polynomial in
    Z_p[x0, ..., x{k-2}, y] ~ Z_p[y][x0, ..., x{k-2}].
    """
    ring = f.ring
    dom = ring.domain
    k = ring.ngens

    coeffs = {}
    for monom, coeff in f.iteritems():
        if not coeffs.has_key(monom[:-1]):
            coeffs[monom[:-1]] = {}
        coeffs[monom[:-1]][monom[-1]] = coeff

    cont = []
    for coeff in coeffs.itervalues():
        cont = gf_gcd(cont, gf_from_dict(coeff, p, dom), p, dom)

    yring = ring.clone(symbols=ring.symbols[k-1])
    contf = yring.from_dense(cont).trunc_ground(p)

    return contf, f.quo(contf.set_ring(ring))


def _LC(f):
    """
    Compute the leading coefficient of a polynomial in
    Z_p[x0, ..., x{k-2}, y] ~ Z_p[y][x0, ..., x{k-2}].
    """
    ring = f.ring
    k = ring.ngens
    yring = ring.clone(symbols=ring.symbols[k-1])
    y = yring.gens[0]
    degf = _deg(f)

    lcf = yring.zero
    for monom, coeff in f.iteritems():
        if monom[:-1] == degf:
            lcf += coeff*y**monom[-1]
    return lcf


def _deg(f):
    """
    Compute the degree of a polynomial in
    Z_p[x0, ..., x{k-2}, y] ~ Z_p[y][x0, ..., x{k-2}].
    """
    k = f.ring.ngens
    degf = tuple(0 for i in range(k-1))
    for monom in f.iterkeys():
        if monom[:-1] > degf:
            degf = monom[:-1]
    return degf


def modgcd_multivariate(f, g):
    """
    Compute the GCD of two polynomials in Z[x0, ..., x{k-1}].
    """
    assert f.ring == g.ring and f.ring.domain.is_ZZ

    result = _trivial_gcd(f, g)
    if result is not None:
        return result

    ring = f.ring
    k = ring.ngens

    # divide out integer content
    cf, f = f.primitive()
    cg, g = g.primitive()
    ch = ring.domain.gcd(cf, cg)

    gamma = ring.domain.gcd(f.LC, g.LC)

    # TODO: write bound function
    degbound = []
    contbound = []
    badprimes = ring.domain.one
    for i in range(k):
        bound = min(ring.dmp_degree_in(f, i), ring.dmp_degree_in(f, i))
        degbound.append(bound)
        contbound.append(bound)

        badprimes *= ring.domain.gcd(_swap(f, i).LC, _swap(g, i).LC)

    m = 1
    p = 1

    while True:
        p = nextprime(p)
        while badprimes % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)

        try:
            hp = modgcd_multivariate_p(f, g, p, degbound, contbound)
            # monic GCD of fp, gp in Z_p[X, y]
        except ModularGCDFailed:
            m = 1
            continue

        if hp is None:
            continue

        hp = hp.mul_ground(gamma).trunc_ground(p)
        if m == 1:
            m = p
            hlastm = hp
            continue

        hm = _chinese_remainder_reconstruction_multivariate(hp, hlastm, p, m)
        m *= p

        if not hm == hlastm:
            hlastm = hm
            continue

        h = hm.primitive()[1]
        fquo, frem = f.div(h)
        gquo, grem = g.div(h)
        if not frem and not grem:
            if h.LC < 0:
                ch = -ch
            h = h.mul_ground(ch)
            cff = fquo.mul_ground(cf // ch)
            cfg = gquo.mul_ground(cg // ch)
            return h, cff, cfg


def modgcd_multivariate_p(f, g, p, degbound, contbound):
    """
    Compute the GCD and cofactors of two polynomials in
    Z_p[x0, ..., x{k-1}].
    """
    ring = f.ring
    k = ring.ngens

    if k == 1:
        h = _gf_gcd(f, g, p).trunc_ground(p)
        (degh,) = h.LM # TODO: use h.degree() instead

        if degh > degbound[0]:
            return None
        if degh < degbound[0]:
            degbound[0] = degh
            raise ModularGCDFailed('all previous primes were unlucky')

        return h

    contf, f = _primitive(f, p)
    contg, g = _primitive(g, p)

    conth = _gf_gcd(contf, contg, p) # polynomial in Z_p[y]
    degyf = ring.dmp_degree_in(f, k-1)
    degyg = ring.dmp_degree_in(g, k-1)
    (degcontf,) = contf.LM # TODO: use contf.degree() instead
    (degcontg,) = contg.LM
    (degconth,) = conth.LM

    if degconth > contbound[k-1]:
        return None
    if degconth < contbound[k-1]:
        contbound[k-1] = degconth
        raise ModularGCDFailed

    lcf = _LC(f)
    lcg = _LC(g)
    delta = _gf_gcd(lcf, lcg, p) # polynomial in Z_p[y]
    (degdelta,) = delta.LM
    y = delta.ring.gens[0]

    zero = tuple(0 for i in range(k-1))
    N = min(degyf - degcontf, degyg - degcontg,
            degbound[k-1] - contbound[k-1] + degdelta) + 1

    if p < N:
        raise ModularGCDFailed

    n = 0
    evalpoints = []
    heval = []

    for a in range(p):
        if not lcf.evaluate(y, a) * lcg.evaluate(y, a) % p:
            continue

        deltaa = delta.evaluate(y, a) % p
        y_ = ring.gens[k-1] # problem: y != y_
        fa = f.evaluate(y_, a).mul_ground(deltaa).trunc_ground(p) # polynomial in Z_p[X]
        ga = g.evaluate(y_, a).mul_ground(deltaa).trunc_ground(p)

        ha = modgcd_multivariate_p(fa, ga, p, degbound, contbound) # polynomial in Z_p[X]

        if ha is None:
            continue

        if ha.LM == zero:
            h = conth.set_ring(ring).trunc_ground(p)
            return h

        evalpoints.append(a)
        heval.append(ha)
        n += 1

        if n == N:
            h = _interpolate_multivariate(evalpoints, heval, ring, p)

            h = _primitive(h, p)[1] * conth.set_ring(ring)
            degyh = ring.dmp_degree_in(h, k-1)

            if degyh > degbound[k-1]:
                return None
            if degyh < degbound[k-1]:
                degbound[k-1] = degyh
                raise ModularGCDFailed

            return h

    return None
