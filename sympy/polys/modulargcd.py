from sympy.ntheory import nextprime
from sympy.ntheory.modular import crt

from sympy.polys.galoistools import (
    gf_gcd, gf_from_dict, gf_gcdex, gf_div, gf_lcm, gf_rem)
from sympy.polys.polyerrors import ModularGCDFailed
from sympy.polys.domains import PolynomialRing

from sympy.core.compatibility import xrange
from sympy.mpmath import sqrt
from sympy import Dummy
import random


def _trivial_gcd(f, g):
    """
    Compute the GCD of two polynomials in trivial cases, i.e. when one
    or both polynomials are zero.
    """
    ring = f.ring

    if not (f or g):
        return ring.zero, ring.zero, ring.zero
    elif not f:
        if g.LC < ring.domain.zero:
            return -g, ring.zero, -ring.one
        else:
            return g, ring.zero, ring.one
    elif not g:
        if f.LC < ring.domain.zero:
            return -f, -ring.one, ring.zero
        else:
            return f, ring.one, ring.zero
    return None


def _gf_gcd(fp, gp, p):
    r"""
    Compute the GCD of two univariate polynomials in `\mathbb{Z}_p[x]`.
    """
    dom = fp.ring.domain

    while gp:
        rem = fp
        deg = gp.degree()
        lcinv = dom.invert(gp.LC, p)

        while True:
            degrem = rem.degree()
            if degrem < deg:
                break
            rem = (rem - gp.mul_monom((degrem - deg,)).mul_ground(lcinv * rem.LC)).trunc_ground(p)

        fp = gp
        gp = rem

    return fp.mul_ground(dom.invert(fp.LC, p)).trunc_ground(p)


def _degree_bound_univariate(f, g):
    r"""
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
    r"""
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
    r"""
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
        cofactor of `f`, i.e. `\frac{f}{h}`
    cfg : PolyElement
        cofactor of `g`, i.e. `\frac{g}{h}`

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


def _primitive(f, p):
    r"""
    Compute the content and the primitive part of a polynomial in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-2}, y] \cong \mathbb{Z}_p[y][x_0, \ldots, x_{k-2}]`.

    Parameters
    ==========

    f : PolyElement
        integer polynomial in `\mathbb{Z}_p[x0, \ldots, x{k-2}, y]`
    p : Integer
        modulus of `f`

    Returns
    =======

    contf : PolyElement
        integer polynomial in `\mathbb{Z}_p[y]`, content of `f`
    ppf : PolyElement
        primitive part of `f`, i.e. `\frac{f}{contf}`

    Examples
    ========

    >>> from sympy.polys.modulargcd import _primitive
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)
    >>> p = 3

    >>> f = x**2*y**2 + x**2*y - y**2 - y
    >>> _primitive(f, p)
    (y**2 + y, x**2 - 1)

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x*y*z - y**2*z**2
    >>> _primitive(f, p)
    (z, x*y - y**2*z)

    """
    ring = f.ring
    dom = ring.domain
    k = ring.ngens

    coeffs = {}
    for monom, coeff in f.iterterms():
        if monom[:-1] not in coeffs:
            coeffs[monom[:-1]] = {}
        coeffs[monom[:-1]][monom[-1]] = coeff

    cont = []
    for coeff in iter(coeffs.values()):
        cont = gf_gcd(cont, gf_from_dict(coeff, p, dom), p, dom)

    yring = ring.clone(symbols=ring.symbols[k-1])
    contf = yring.from_dense(cont).trunc_ground(p)

    return contf, f.quo(contf.set_ring(ring))


def _deg(f):
    r"""
    Compute the degree of a multivariate polynomial
    `f \in K[x_0, \ldots, x_{k-2}, y] \cong K[y][x_0, \ldots, x_{k-2}]`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `K[x_0, \ldots, x_{k-2}, y]`

    Returns
    =======

    degf : Integer tuple
        degree of `f` in `x_0, \ldots, x_{k-2}`

    Examples
    ========

    >>> from sympy.polys.modulargcd import _deg
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _deg(f)
    (2,)

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _deg(f)
    (2, 2)

    >>> f = x*y*z - y**2*z**2
    >>> _deg(f)
    (1, 1)

    """
    k = f.ring.ngens
    degf = (0,) * (k-1)
    for monom in f.itermonoms():
        if monom[:-1] > degf:
            degf = monom[:-1]
    return degf


def _LC(f):
    r"""
    Compute the leading coefficient of a multivariate polynomial
    `f \in K[x_0, \ldots, x_{k-2}, y] \cong K[y][x_0, \ldots, x_{k-2}]`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `K[x_0, \ldots, x_{k-2}, y]`

    Returns
    =======

    lcf : PolyElement
        polynomial in `K[y]`, leading coefficient of `f`

    Examples
    ========

    >>> from sympy.polys.modulargcd import _LC
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _LC(f)
    y**2 + y

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x**2*y**2 + x**2*y - 1
    >>> _LC(f)
    1

    >>> f = x*y*z - y**2*z**2
    >>> _LC(f)
    z

    """
    ring = f.ring
    k = ring.ngens
    yring = ring.clone(symbols=ring.symbols[k-1])
    y = yring.gens[0]
    degf = _deg(f)

    lcf = yring.zero
    for monom, coeff in f.iterterms():
        if monom[:-1] == degf:
            lcf += coeff*y**monom[-1]
    return lcf


def _swap(f, i):
    """
    Make the variable `x_i` the leading one in a multivariate polynomial `f`.
    """
    ring = f.ring
    k = ring.ngens
    fswap = ring.zero
    for monom, coeff in f.iterterms():
        monomswap = (monom[i],) + monom[:i] + monom[i+1:]
        fswap[monomswap] = coeff
    return fswap


def _degree_bound_bivariate(f, g):
    r"""
    Compute upper degree bounds for the GCD of two bivariate
    integer polynomials `f` and `g`.

    The GCD is viewed as a polynomial in `\mathbb{Z}[y][x]` and the
    function returns an upper bound for its degree and one for the degree
    of its content. This is done by choosing a suitable prime `p` and
    computing the GCD of the contents of `f \; \mathrm{mod} \, p` and
    `g \; \mathrm{mod} \, p`. The choice of `p` guarantees that the degree
    of the content in `\mathbb{Z}_p[y]` is greater than or equal to the
    degree in `\mathbb{Z}[y]`. To obtain the degree bound in the variable
    `x`, the polynomials are evaluated at `y = a` for a suitable
    `a \in \mathbb{Z}_p` and then their GCD in `\mathbb{Z}_p[x]` is
    computed. If no such `a` exists, i.e. the degree in `\mathbb{Z}_p[x]`
    is always smaller than the one in `\mathbb{Z}[y][x]`, then the bound is
    set to the minimum of the degrees of `f` and `g` in `x`.

    Parameters
    ==========

    f : PolyElement
        bivariate integer polynomial
    g : PolyElement
        bivariate integer polynomial

    Returns
    =======

    xbound : Integer
        upper bound for the degree of the GCD of the polynomials `f` and
        `g` in the variable `x`
    ycontbound : Integer
        upper bound for the degree of the content of the GCD of the
        polynomials `f` and `g` in the variable `y`

    References
    ==========

    1. [Monagan00]_

    """
    ring = f.ring

    gamma1 = ring.domain.gcd(f.LC, g.LC)
    gamma2 = ring.domain.gcd(_swap(f, 1).LC, _swap(g, 1).LC)
    badprimes = gamma1 * gamma2
    p = 1

    while True:
        p = nextprime(p)
        while badprimes % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        contfp, fp = _primitive(fp, p)
        contgp, gp = _primitive(gp, p)
        conthp = _gf_gcd(contfp, contgp, p) # polynomial in Z_p[y]
        ycontbound = conthp.degree()

        # polynomial in Z_p[y]
        delta = _gf_gcd(_LC(fp), _LC(gp), p)

        for a in xrange(p):
            if not delta.evaluate(0, a) % p:
                continue
            fpa = fp.evaluate(1, a).trunc_ground(p)
            gpa = gp.evaluate(1, a).trunc_ground(p)
            hpa = _gf_gcd(fpa, gpa, p)
            xbound = hpa.degree()
            return xbound, ycontbound

        return min(fp.degree(), gp.degree()), ycontbound


def _chinese_remainder_reconstruction_multivariate(hp, hq, p, q):
    r"""
    Construct a polynomial `h_{pq}` in
    `\mathbb{Z}_{p q}[x_0, \ldots, x_{k-1}]` such that

    .. math ::

        h_{pq} = h_p \; \mathrm{mod} \, p

        h_{pq} = h_q \; \mathrm{mod} \, q

    for relatively prime integers `p` and `q` and polynomials
    `h_p` and `h_q` in `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]` and
    `\mathbb{Z}_q[x_0, \ldots, x_{k-1}]` respectively.

    The coefficients of the polynomial `h_{pq}` are computed with the
    Chinese Remainder Theorem. The symmetric representation in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`,
    `\mathbb{Z}_q[x_0, \ldots, x_{k-1}]` and
    `\mathbb{Z}_{p q}[x_0, \ldots, x_{k-1}]` is used.

    Parameters
    ==========

    hp : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    hq : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_q`
    p : Integer
        modulus of `h_p`, relatively prime to `q`
    q : Integer
        modulus of `h_q`, relatively prime to `p`

    Examples
    ========

    >>> from sympy.polys.modulargcd import _chinese_remainder_reconstruction_multivariate
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)
    >>> p = 3
    >>> q = 5

    >>> hp = x**3*y - x**2 - 1
    >>> hq = -x**3*y - 2*x*y**2 + 2

    >>> hpq = _chinese_remainder_reconstruction_multivariate(hp, hq, p, q)
    >>> hpq
    4*x**3*y + 5*x**2 + 3*x*y**2 + 2

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    >>> R, x, y, z = ring("x, y, z", ZZ)
    >>> p = 6
    >>> q = 5

    >>> hp = 3*x**4 - y**3*z + z
    >>> hq = -2*x**4 + z

    >>> hpq = _chinese_remainder_reconstruction_multivariate(hp, hq, p, q)
    >>> hpq
    3*x**4 + 5*y**3*z + z

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    """
    hpmonoms = set(hp.monoms())
    hqmonoms = set(hq.monoms())
    monoms = hpmonoms.intersection(hqmonoms)
    hpmonoms.difference_update(monoms)
    hqmonoms.difference_update(monoms)

    zero = hp.ring.domain.zero

    hpq = hp.ring.zero

    if isinstance(hp.ring.domain, PolynomialRing):
        crt_ = _chinese_remainder_reconstruction_multivariate
    else:
        def crt_(cp, cq, p, q):
            return crt([p, q], [cp, cq], symmetric=True)[0]

    for monom in monoms:
        hpq[monom] = crt_(hp[monom], hq[monom], p, q)
    for monom in hpmonoms:
        hpq[monom] = crt_(hp[monom], zero, p, q)
    for monom in hqmonoms:
        hpq[monom] = crt_(zero, hq[monom], p, q)

    return hpq


def _interpolate_multivariate(evalpoints, hpeval, ring, i, p, ground=False):
    r"""
    Reconstruct a polynomial `h_p` in `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`
    from a list of evaluation points in `\mathbb{Z}_p` and a list of
    polynomials in `\mathbb{Z}_p[x_0, \ldots, x_{k-2}]`, which are the images
    of `h_p` evaluated in the variable `x_{k-1}`.

    Parameters
    ==========

    evalpoints : list of Integer objects
        list of evaluation points in `\mathbb{Z}_p`
    hpeval : list of PolyElement objects
        list of polynomials in `\mathbb{Z}_p[x_0, \ldots, x_{k-2}]`,
        images of `h_p` evaluated in the variable `x_{k-1}`
    ring : PolyRing
        `h_p` will be an element of this ring
    p : Integer
        prime number, modulus of `h_p`

    Returns
    =======

    hp : PolyElement
        interpolated polynomial in `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`

    """
    hp = ring.zero

    if ground:
        domain = ring.domain.domain
        y = ring.domain.gens[i]
    else:
        domain = ring.domain
        y = ring.gens[i]

    for a, hpa in zip(evalpoints, hpeval):
        numer = ring.one
        denom = domain.one
        for b in evalpoints:
            if b == a:
                continue

            numer *= y - b
            denom *= a - b

        denom = domain.invert(denom, p)
        coeff = numer.mul_ground(denom)
        hp += hpa.set_ring(ring) * coeff

    return hp.trunc_ground(p)


def modgcd_bivariate(f, g):
    r"""
    Computes the GCD of two polynomials in `\mathbb{Z}[x, y]` using a
    modular algorithm.

    The algorithm computes the GCD of two bivariate integer polynomials
    `f` and `g` by calculating the GCD in `\mathbb{Z}_p[x, y]` for
    suitable primes `p` and then reconstructing the coefficients with the
    Chinese Remainder Theorem. To compute the bivariate GCD over
    `\mathbb{Z}_p`, the polynomials `f \; \mathrm{mod} \, p` and
    `g \; \mathrm{mod} \, p` are evaluated at `y = a` for certain
    `a \in \mathbb{Z}_p` and then their univariate GCD in `\mathbb{Z}_p[x]`
    is computed. Interpolating those yields the bivariate GCD in
    `\mathbb{Z}_p[x, y]`. To verify the result in `\mathbb{Z}[x, y]`, trial
    division is done, but only for candidates which are very likely the
    desired GCD.

    Parameters
    ==========

    f : PolyElement
        bivariate integer polynomial
    g : PolyElement
        bivariate integer polynomial

    Returns
    =======

    h : PolyElement
        GCD of the polynomials `f` and `g`
    cff : PolyElement
        cofactor of `f`, i.e. `\frac{f}{h}`
    cfg : PolyElement
        cofactor of `g`, i.e. `\frac{g}{h}`

    Examples
    ========

    >>> from sympy.polys.modulargcd import modgcd_bivariate
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2 - y**2
    >>> g = x**2 + 2*x*y + y**2

    >>> h, cff, cfg = modgcd_bivariate(f, g)
    >>> h, cff, cfg
    (x + y, x - y, x + y)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    >>> f = x**2*y - x**2 - 4*y + 4
    >>> g = x + 2

    >>> h, cff, cfg = modgcd_bivariate(f, g)
    >>> h, cff, cfg
    (x + 2, x*y - x - 2*y + 2, 1)

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

    xbound, ycontbound = _degree_bound_bivariate(f, g)
    if xbound == ycontbound == 0:
        return ring(ch), f.mul_ground(cf // ch), g.mul_ground(cg // ch)

    fswap = _swap(f, 1)
    gswap = _swap(g, 1)
    degyf = fswap.degree()
    degyg = gswap.degree()

    ybound, xcontbound = _degree_bound_bivariate(fswap, gswap)
    if ybound == xcontbound == 0:
        return ring(ch), f.mul_ground(cf // ch), g.mul_ground(cg // ch)

    # TODO: to improve performance, choose the main variable here

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
        contfp, fp = _primitive(fp, p)
        contgp, gp = _primitive(gp, p)
        conthp = _gf_gcd(contfp, contgp, p) # monic polynomial in Z_p[y]
        degconthp = conthp.degree()

        if degconthp > ycontbound:
            continue
        elif degconthp < ycontbound:
            m = 1
            ycontbound = degconthp
            continue

        # polynomial in Z_p[y]
        delta = _gf_gcd(_LC(fp), _LC(gp), p)

        degcontfp = contfp.degree()
        degcontgp = contgp.degree()
        degdelta = delta.degree()

        N = min(degyf - degcontfp, degyg - degcontgp,
            ybound - ycontbound + degdelta) + 1

        if p < N:
            continue

        n = 0
        evalpoints = []
        hpeval = []
        unlucky = False

        for a in xrange(p):
            deltaa = delta.evaluate(0, a)
            if not deltaa % p:
                continue

            fpa = fp.evaluate(1, a).trunc_ground(p)
            gpa = gp.evaluate(1, a).trunc_ground(p)
            hpa = _gf_gcd(fpa, gpa, p) # monic polynomial in Z_p[x]
            deghpa = hpa.degree()

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

        hp = _interpolate_multivariate(evalpoints, hpeval, ring, 1, p)

        hp = _primitive(hp, p)[1]
        hp = hp * conthp.set_ring(ring)
        degyhp = hp.degree(1)

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


def _modgcd_multivariate_p(f, g, p, degbound, contbound):
    r"""
    Compute the GCD of two polynomials in
    `\mathbb{Z}_p[x0, \ldots, x{k-1}]`.

    The algorithm reduces the problem step by step by evaluating the
    polynomials `f` and `g` at `x_{k-1} = a` for suitable
    `a \in \mathbb{Z}_p` and then calls itself recursively to compute the GCD
    in `\mathbb{Z}_p[x_0, \ldots, x_{k-2}]`. If these recursive calls are
    succsessful for enough evaluation points, the GCD in `k` variables is
    interpolated, otherwise the algorithm returns ``None``. Every time a GCD
    or a content is computed, their degrees are compared with the bounds. If
    a degree greater then the bound is encountered, then the current call
    returns ``None`` and a new evaluation point has to be chosen. If at some
    point the degree is smaller, the correspondent bound is updated and the
    algorithm fails.

    Parameters
    ==========

    f : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    g : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    p : Integer
        prime number, modulus of `f` and `g`
    degbound : list of Integer objects
        ``degbound[i]`` is an upper bound for the degree of the GCD of `f`
        and `g` in the variable `x_i`
    contbound : list of Integer objects
        ``contbound[i]`` is an upper bound for the degree of the content of
        the GCD in `\mathbb{Z}_p[x_i][x_0, \ldots, x_{i-1}]`,
        ``contbound[0]`` is not used can therefore be chosen
        arbitrarily.

    Returns
    =======

    h : PolyElement
        GCD of the polynomials `f` and `g` or ``None``

    References
    ==========

    1. [Monagan00]_
    2. [Brown71]_

    """
    ring = f.ring
    k = ring.ngens

    if k == 1:
        h = _gf_gcd(f, g, p).trunc_ground(p)
        degh = h.degree()

        if degh > degbound[0]:
            return None
        if degh < degbound[0]:
            degbound[0] = degh
            raise ModularGCDFailed

        return h

    degyf = f.degree(k-1)
    degyg = g.degree(k-1)

    contf, f = _primitive(f, p)
    contg, g = _primitive(g, p)

    conth = _gf_gcd(contf, contg, p) # polynomial in Z_p[y]

    degcontf = contf.degree()
    degcontg = contg.degree()
    degconth = conth.degree()

    if degconth > contbound[k-1]:
        return None
    if degconth < contbound[k-1]:
        contbound[k-1] = degconth
        raise ModularGCDFailed

    lcf = _LC(f)
    lcg = _LC(g)

    delta = _gf_gcd(lcf, lcg, p) # polynomial in Z_p[y]

    evaltest = delta

    for i in xrange(k-1):
        evaltest *= _gf_gcd(_LC(_swap(f, i)), _LC(_swap(g, i)), p)

    degdelta = delta.degree()

    N = min(degyf - degcontf, degyg - degcontg,
            degbound[k-1] - contbound[k-1] + degdelta) + 1

    if p < N:
        return None

    n = 0
    d = 0
    evalpoints = []
    heval = []
    points = set(range(p))

    while points:
        a = random.sample(points, 1)[0]
        points.remove(a)

        if not evaltest.evaluate(0, a) % p:
            continue

        deltaa = delta.evaluate(0, a) % p

        fa = f.evaluate(k-1, a).trunc_ground(p)
        ga = g.evaluate(k-1, a).trunc_ground(p)

        # polynomials in Z_p[x_0, ..., x_{k-2}]
        ha = _modgcd_multivariate_p(fa, ga, p, degbound, contbound)

        if ha is None:
            d += 1
            if d > n:
                return None
            continue

        if ha.is_ground:
            h = conth.set_ring(ring).trunc_ground(p)
            return h

        ha = ha.mul_ground(deltaa).trunc_ground(p)

        evalpoints.append(a)
        heval.append(ha)
        n += 1

        if n == N:
            h = _interpolate_multivariate(evalpoints, heval, ring, k-1, p)

            h = _primitive(h, p)[1] * conth.set_ring(ring)
            degyh = h.degree(k-1)

            if degyh > degbound[k-1]:
                return None
            if degyh < degbound[k-1]:
                degbound[k-1] = degyh
                raise ModularGCDFailed

            return h

    return None


def modgcd_multivariate(f, g):
    r"""
    Compute the GCD of two polynomials in `\mathbb{Z}[x_0, \ldots, x_{k-1}]`
    using a modular algorithm.

    The algorithm computes the GCD of two multivariate integer polynomials
    `f` and `g` by calculating the GCD in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]` for suitable primes `p` and then
    reconstructing the coefficients with the Chinese Remainder Theorem. To
    compute the multivariate GCD over `\mathbb{Z}_p` the recursive
    subroutine ``_modgcd_multivariate_p`` is used. To verify the result in
    `\mathbb{Z}[x_0, \ldots, x_{k-1}]`, trial division is done, but only for
    candidates which are very likely the desired GCD.

    Parameters
    ==========

    f : PolyElement
        multivariate integer polynomial
    g : PolyElement
        multivariate integer polynomial

    Returns
    =======

    h : PolyElement
        GCD of the polynomials `f` and `g`
    cff : PolyElement
        cofactor of `f`, i.e. `\frac{f}{h}`
    cfg : PolyElement
        cofactor of `g`, i.e. `\frac{g}{h}`

    Examples
    ========

    >>> from sympy.polys.modulargcd import modgcd_multivariate
    >>> from sympy.polys import ring, ZZ

    >>> R, x, y = ring("x, y", ZZ)

    >>> f = x**2 - y**2
    >>> g = x**2 + 2*x*y + y**2

    >>> h, cff, cfg = modgcd_multivariate(f, g)
    >>> h, cff, cfg
    (x + y, x - y, x + y)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    >>> R, x, y, z = ring("x, y, z", ZZ)

    >>> f = x*z**2 - y*z**2
    >>> g = x**2*z + z

    >>> h, cff, cfg = modgcd_multivariate(f, g)
    >>> h, cff, cfg
    (z, x*z - y*z, x**2 + 1)

    >>> cff * h == f
    True
    >>> cfg * h == g
    True

    References
    ==========

    1. [Monagan00]_
    2. [Brown71]_

    See also
    ========

    _modgcd_multivariate_p

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

    badprimes = ring.domain.one
    for i in xrange(k):
        badprimes *= ring.domain.gcd(_swap(f, i).LC, _swap(g, i).LC)

    degbound = [min(fdeg, gdeg) for fdeg, gdeg in zip(f.degrees(), g.degrees())]
    contbound = list(degbound)

    m = 1
    p = 1

    while True:
        p = nextprime(p)
        while badprimes % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)

        try:
            # monic GCD of fp, gp in Z_p[x_0, ..., x_{k-2}, y]
            hp = _modgcd_multivariate_p(fp, gp, p, degbound, contbound)
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


def _gf_div(f, g, p):
    ring = f.ring
    densequo, denserem = gf_div(f.to_dense(), g.to_dense(), p, ring.domain)
    return ring.from_dense(densequo), ring.from_dense(denserem)


def _rational_function_reconstruction(c, p, m):
    ring = c.ring
    domain = ring.domain
    M = m.degree()
    N = M // 2
    D = M - N - 1

    r0, s0 = m, ring.zero
    r1, s1 = c, ring.one

    while r1.degree() > N:
        quo = _gf_div(r0, r1, p)[0]
        r0, r1 = r1, (r0 - quo*r1).trunc_ground(p)
        s0, s1 = s1, (s0 - quo*s1).trunc_ground(p)

    a, b = r1, s1
    if b.degree() > D or _gf_gcd(b, m, p) != 1:
        return None

    lc = b.LC
    if lc != 1:
        lcinv = domain.invert(lc, p)
        a = a.mul_ground(lcinv).trunc_ground(p)
        b = b.mul_ground(lcinv).trunc_ground(p)

    field = ring.to_field()

    return field(a) / field(b)


def _rational_reconstruction_func_coeffs(hm, p, m, ring, k):
    h = ring.zero

    for monom, coeff in hm.iterterms():
        if k == 0:
            coeffh = _rational_function_reconstruction(coeff, p, m)

            if not coeffh:
                return None

        else:
            coeffh = ring.domain.zero
            for mon, c in coeff.drop_to_ground(k).iterterms():
                ch = _rational_function_reconstruction(c, p, m)

                if not ch:
                    return None

                coeffh[mon] = ch

        h[monom] = coeffh

    return h


def _gf_gcdex(f, g, p):
    ring = f.ring
    s, t, h = gf_gcdex(f.to_dense(), g.to_dense(), p, ring.domain)
    return ring.from_dense(s), ring.from_dense(t), ring.from_dense(h)


def _trunc(f, minpoly, p):



    ring = f.ring
    minpoly = minpoly.set_ring(ring)
    p_ = ring.ground_new(p)

    return f.trunc_ground(p).rem([minpoly, p_]).trunc_ground(p)


def _euclidean_algorithm(f, g, minpoly, p):
    ring = f.ring

    f = _trunc(f, minpoly, p)
    g = _trunc(g, minpoly, p)

    while g:
        rem = f
        deg = g.degree(0) # degree in x
        lcinv, _, gcd = _gf_gcdex(ring.dmp_LC(g), minpoly, p)

        if not gcd == 1:
            return None

        while True:
            degrem = rem.degree(0) # degree in x
            if degrem < deg:
                break
            quo = (lcinv * ring.dmp_LC(rem)).set_ring(ring)
            rem = _trunc(rem - g.mul_monom((degrem - deg, 0))*quo, minpoly, p)

        f = g
        g = rem

    lcfinv = _gf_gcdex(ring.dmp_LC(f), minpoly, p)[0].set_ring(ring)

    return _trunc(f * lcfinv, minpoly, p)


def _trial_division(f, h, minpoly, p=None):
    ring = f.ring
    domain = ring.domain

    zxring = ring.clone(symbols=(ring.symbols[1], ring.symbols[0]))

    minpoly = minpoly.set_ring(ring)

    rem = f

    degrem = rem.degree()
    degh = h.degree()
    degm = minpoly.degree(1)

    lch = _LC(h).set_ring(ring)
    lcm = minpoly.LC

    while rem and degrem >= degh:
        # polynomial in Z[t_1, ..., t_k][z]
        lcrem = _LC(rem).set_ring(ring)
        rem = rem*lch - h.mul_monom((degrem - degh, 0))*lcrem
        if p:
            rem = rem.trunc_ground(p)
        degrem = rem.degree(1)

        while rem and degrem >= degm:
            # polynomial in Z[t_1, ..., t_k][x]
            lcrem = _LC(rem.set_ring(zxring)).set_ring(ring)
            rem = rem.mul_ground(lcm) - minpoly.mul_monom((0, degrem - degm))*lcrem
            if p:
                rem = rem.trunc_ground(p)
            degrem = rem.degree(1)

        degrem = rem.degree()

    return rem


def _evaluate_ground(f, i, a):
    ring = f.ring.clone(domain=f.ring.domain.ring.drop(i))
    fa = ring.zero

    for monom, coeff in f.iterterms():
        fa[monom] = coeff.evaluate(i, a)

    return fa


def _func_field_modgcd_p(f, g, minpoly, p):
    ring = f.ring
    domain = ring.domain # Z[t_1, ..., t_k]

    if isinstance(domain, PolynomialRing):
        k = domain.ngens
    else:
        return _euclidean_algorithm(f, g, minpoly, p)

    if k == 1:
        qdomain = domain.ring.to_field()
    else:
        qdomain = domain.ring.drop_to_ground(k - 1)
        qdomain = qdomain.clone(domain=qdomain.domain.ring.to_field())

    qring = ring.clone(domain=qdomain) # = Z(t_k)[t_1, ..., t_{k-1}][x, z]

    n = 1
    d = 1

    # polynomial in Z_p[t_1, ..., t_k][z]
    gamma = ring.dmp_LC(f) * ring.dmp_LC(g)
    # polynomial in Z_p[t_1, ..., t_k]
    delta = minpoly.LC

    evalpoints = []
    heval = []
    LMlist = []
    points = set(range(p))

    while points:
        a = random.sample(points, 1)[0]
        points.remove(a)

        if k == 1:
            test = delta.evaluate(k-1, a) % p == 0
        else:
            test = delta.evaluate(k-1, a).trunc_ground(p) == 0

        if test:
            continue

        gammaa = _evaluate_ground(gamma, k-1, a)
        minpolya = _evaluate_ground(minpoly, k-1, a)

        if gammaa.rem([minpolya, gammaa.ring(p)]) == 0:
            continue

        fa = _evaluate_ground(f, k-1, a)
        ga = _evaluate_ground(g, k-1, a)

        # polynomial in Z_p[x, t_1, ..., t_{k-1}, z]/(minpoly)
        ha = _func_field_modgcd_p(fa, ga, minpolya, p)

        if ha is None:
            d += 1
            if d > n:
                return None
            continue

        if ha == 1:
            return ha

        LM = [ha.degree()] + [0]*(k-1)
        if k > 1:
            for monom, coeff in ha.iterterms():
                if monom[0] == LM[0] and coeff.LM > tuple(LM[1:]):
                    LM[1:] = coeff.LM

        evalpoints_a = [a]
        heval_a = [ha]
        if k == 1:
            m = qring.domain.get_ring().one
        else:
            m = qring.domain.domain.get_ring().one

        t = m.ring.gens[0]

        for b, hb, LMhb in zip(evalpoints, heval, LMlist):
            if LMhb == LM:
                evalpoints_a.append(b)
                heval_a.append(hb)
                m *= (t - b)

        m = m.trunc_ground(p)
        evalpoints.append(a)
        heval.append(ha)
        LMlist.append(LM)
        n += 1

        # polynomial in Z_p[t_1, ..., t_k][x, z]
        h = _interpolate_multivariate(evalpoints_a, heval_a, ring, k-1, p, ground=True)

        # polynomial in Z_p(t_k)[t_1, ..., t_{k-1}][x, z]
        h = _rational_reconstruction_func_coeffs(h, p, m, qring, k-1)

        if h is None:
            continue

        if k == 1:
            dom = qring.domain.field
            den = dom.ring.one

            for coeff in h.itercoeffs():
                den = dom.ring.from_dense(gf_lcm(den.to_dense(), coeff.denom.to_dense(),
                        p, dom.domain))

        else:
            dom = qring.domain.domain.field
            den = dom.ring.one

            for coeff in h.itercoeffs():
                for c in coeff.itercoeffs():
                    den = dom.ring.from_dense(gf_lcm(den.to_dense(), c.denom.to_dense(),
                            p, dom.domain))

        den = qring.domain_new(den.trunc_ground(p))
        h = ring(h.mul_ground(den).as_expr()).trunc_ground(p)

        if not _trial_division(f, h, minpoly, p) and not _trial_division(g, h, minpoly, p):
            return h

    return None


def _integer_rational_reconstruction(c, m, domain):
    if c < 0:
        c += m

    r0, s0 = m, domain.zero
    r1, s1 = c, domain.one

    bound = sqrt(m / 2) # still correct if replaced by ZZ.sqrt(m // 2) ?

    while r1 >= bound:
        quo = r0 // r1
        r0, r1 = r1, r0 - quo*r1
        s0, s1 = s1, s0 - quo*s1

    if abs(s1) >= bound:
        return None

    if s1 < 0:
        a, b = -r1, -s1
    elif s1 > 0:
        a, b = r1, s1
    else:
        return None

    field = domain.get_field()

    return field(a) / field(b)


def _rational_reconstruction_int_coeffs(hm, m, ring):
    h = ring.zero

    if isinstance(ring.domain, PolynomialRing):
        reconstruction = _rational_reconstruction_int_coeffs
        domain = ring.domain.ring
    else:
        reconstruction = _integer_rational_reconstruction
        domain = hm.ring.domain

    for monom, coeff in hm.iterterms():
        coeffh = reconstruction(coeff, m, domain)

        if not coeffh:
            return None

        h[monom] = coeffh

    return h


def _func_field_modgcd_m(f, g, minpoly):
    ring = f.ring
    domain = ring.domain

    if isinstance(domain, PolynomialRing):
        k = domain.ngens
        QQdomain = domain.ring.clone(domain=domain.domain.get_field())
        QQring = ring.clone(domain=QQdomain)
    else:
        k = 0
        QQring = ring.clone(domain=ring.domain.get_field())

    cf, f = f.primitive()
    cg, g = g.primitive()

    # polynomial in Z[t_1, ..., t_k][z]
    gamma = ring.dmp_LC(f) * ring.dmp_LC(g)
    # polynomial in Z[t_1, ..., t_k]
    delta = minpoly.LC

    p = 1
    primes = []
    hplist = []
    LMlist = []

    while True:
        p = nextprime(p)

        if gamma.trunc_ground(p) == 0:
            continue

        if k == 0:
            test = (delta % p == 0)
        else:
            test = (delta.trunc_ground(p) == 0)

        if test:
            continue

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        minpolyp = minpoly.trunc_ground(p)

        hp = _func_field_modgcd_p(fp, gp, minpolyp, p)

        if hp is None:
            continue

        if hp == 1:
            return ring.one

        LM = [hp.degree()] + [0]*k
        if k > 0:
            for monom, coeff in hp.iterterms():
                if monom[0] == LM[0] and coeff.LM > tuple(LM[1:]):
                    LM[1:] = coeff.LM

        hm = hp
        m = p

        for q, hq, LMhq in zip(primes, hplist, LMlist):
            if LMhq == LM:
                hm = _chinese_remainder_reconstruction_multivariate(hq, hm, q, m)
                m *= q

        primes.append(p)
        hplist.append(hp)
        LMlist.append(LM)

        hm = _rational_reconstruction_int_coeffs(hm, m, QQring)

        if hm is None:
            continue

        if k == 0:
            h = hm.clear_denoms()[1]
        else:
            den = domain.domain.one
            for coeff in hm.itercoeffs():
                den = domain.domain.lcm(den, coeff.clear_denoms()[0])
            h = hm.mul_ground(den)

        # convert back to Z[t_1, ..., t_k][x, z] from Q[t_1, ..., t_k][x, z]
        h = h.set_ring(ring)
        h = h.primitive()[1]

        if not (_trial_division(f.mul_ground(cf), h, minpoly) or
            _trial_division(g.mul_ground(cg), h, minpoly)):
            return h


def _to_ZZ_poly(f, ring):
    f_ = ring.zero

    if isinstance(ring.domain, PolynomialRing):
        domain = ring.domain.domain
    else:
        domain = ring.domain

    den = domain.one

    for coeff in f.itercoeffs():
        for c in coeff.rep:
            if c:
                den = domain.lcm(den, c.denominator)

    for monom, coeff in f.iterterms():
        coeff = coeff.rep
        m = ring.domain.one
        if isinstance(ring.domain, PolynomialRing):
            m = m.mul_monom(monom[1:])
        n = len(coeff)

        for i in xrange(n):
            if coeff[i]:
                c = domain(coeff[i] * den) * m

                if (monom[0], n-i-1) not in f_:
                    f_[(monom[0], n-i-1)] = c
                else:
                    f_[(monom[0], n-i-1)] += c

    return f_


def _to_ANP_poly(f, ring):
    domain = ring.domain
    f_ = ring.zero

    if isinstance(f.ring.domain, PolynomialRing):
        for monom, coeff in f.iterterms():
            for mon, coef in coeff.iterterms():
                m = (monom[0],) + mon
                c = domain([domain.domain(coef)] + [0]*monom[1])

                if m not in f_:
                    f_[m] = c
                else:
                    f_[m] += c

    else:
        for monom, coeff in f.iterterms():
            m = (monom[0],)
            c = domain([domain.domain(coeff)] + [0]*monom[1])

            if m not in f_:
                f_[m] = c
            else:
                f_[m] += c

    return f_


def _minpoly_from_dense(minpoly, ring):
    minpoly_ = ring.zero

    for monom, coeff in minpoly.terms():
        minpoly_[monom] = ring.domain(coeff)

    return minpoly_


def _primitive_in_x0(f):
    fring = f.ring
    ring = fring.drop_to_ground(*xrange(1, fring.ngens))
    dom = ring.domain.ring
    f_ = ring(f.as_expr())
    cont = dom.zero

    for coeff in f_.itercoeffs():
        cont = func_field_modgcd(cont, coeff)[0]
        if cont == dom.one:
            return cont, f

    return cont, f.quo(cont.set_ring(fring))


# TODO: add support for algebraic function fields
def func_field_modgcd(f, g):
    ring = f.ring
    domain = ring.domain
    n = ring.ngens

    assert ring == g.ring and domain.is_Algebraic

    result = _trivial_gcd(f, g)
    if result is not None:
        return result

    z = Dummy('z')

    ZZring = ring.clone(symbols=ring.symbols + (z,), domain=domain.domain.get_ring())

    if n == 1:
        f_ = _to_ZZ_poly(f, ZZring)
        g_ = _to_ZZ_poly(g, ZZring)
        minpoly = ZZring.drop(0).from_dense(domain.mod.rep)

        h = _func_field_modgcd_m(f_, g_, minpoly)
        h = _to_ANP_poly(h, ring)

    else:
        # contx0f in Q(a)[x_1, ..., x_{n-1}], f in Q(a)[x_0, ..., x_{n-1}]
        contx0f, f = _primitive_in_x0(f)
        contx0g, g = _primitive_in_x0(g)
        contx0h = func_field_modgcd(contx0f, contx0g)[0]

        ZZring_ = ZZring.drop_to_ground(*xrange(1, n))

        f_ = _to_ZZ_poly(f, ZZring_)
        g_ = _to_ZZ_poly(g, ZZring_)
        minpoly = _minpoly_from_dense(domain.mod, ZZring_.drop(0))

        h = _func_field_modgcd_m(f_, g_, minpoly)
        h = _to_ANP_poly(h, ring)

        contx0h_, h = _primitive_in_x0(h)
        h *= contx0h.set_ring(ring)
        f *= contx0f.set_ring(ring)
        g *= contx0g.set_ring(ring)

    h = h.quo_ground(h.LC)

    return h, f.quo(h), g.quo(h)
