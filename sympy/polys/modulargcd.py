from sympy.ntheory import nextprime
from sympy.ntheory.modular import crt
from sympy.polys.galoistools import gf_gcd, gf_degree

def modgcd(f, g):
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
    (n,) = hp.LM # TODO: use hp.degree() instead
    x = hp.ring.gens[0]
    hpq = hp.ring.zero

    for i in range(n+1):
        hpq[(i,)] = crt([p, q], [hp.coeff(x**i), hq.coeff(x**i)], symmetric=True)[0]

    hpq.strip_zero()
    return hpq


# should I implement a GCD in GF(p) for PolyElement as well?
def _gf_gcd(fp, gp, p):
    ring = fp.ring
    densehp = gf_gcd(fp.to_dense(), gp.to_dense(), p, ring.domain)
    return ring.from_dense(densehp)


def _trivial_gcd(f, g):
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
