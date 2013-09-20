from __future__ import print_function

from sympy import Dummy
import random

from sympy.ntheory import nextprime
from sympy.polys.galoistools import gf_sqf_p, gf_irreducible_p
from sympy.polys.modulargcd import _trunc, _gf_gcdex, _minpoly_from_dense, _euclidean_algorithm
from sympy.polys.polyclasses import ANP
from sympy.polys.polyutils import _sort_factors
from sympy.polys.rings import PolyRing


# TODO
# ====

# -) efficiency of _factor can be improved for irreducible polynomials if the
#    univariate factorization is done before the LC is factored


def _alpha_to_z(f, ring):
    if isinstance(f, ANP):
        ring = ring.drop(*ring.gens[:-1])
        f_ = ring.from_dense(f.rep)

    else:
        f_ = ring.zero
        domain = ring.domain

        for monom, coeff in f.iterterms():
            coeff = coeff.rep
            n = len(coeff)

            for i in xrange(n):
                m = monom + (n-i-1,)
                if coeff[i]:
                    if m not in f_:
                        f_[m] = coeff[i]
                    else:
                        f_[m] += coeff[i]

    return f_


def _z_to_alpha(f, ring):
    domain = ring.domain

    if f.ring.is_univariate:
        f_ = domain.zero
        for (monom,), coeff in f.iterterms():
            f_ += domain([domain.domain(coeff)] + [0]*monom)

    else:
        f_ = ring.zero
        for monom, coeff in f.iterterms():
            m = monom[:-1]
            c = domain([domain.domain(coeff)] + [0]*monom[-1])

            if m not in f_:
                f_[m] = c
            else:
                f_[m] += c

    return f_


def _distinct_prime_divisors(A, domain):
    gcd = domain.gcd
    divisors = []

    for i, a in enumerate(A):
        divisors.append(a)

        for j in xrange(i-1):
            g = gcd(divisors[i], divisors[j])
            divisors[i] = divisors[i] // g
            divisors[j] = divisors[j] // g
            g1 = gcd(divisors[i], g)
            g2 = gcd(divisors[j], g)

            while g1 != 1:
                g1 = gcd(divisors[i], g1)
                divisors[i] = divisors[i] // g1

            while g2 != 1:
                g2 = gcd(divisors[j], g2)
                divisors[j] = divisors[j] // g2

            if divisors[i] == 1 or divisors[j] == 1:
                return None

    return divisors


def _denominator(f, qring):
    f = _alpha_to_z(f, qring)

    ring = qring.domain.get_ring()
    lcm = ring.lcm
    den = ring.one

    for coeff in f.itercoeffs():
        den = lcm(den, coeff.denominator)

    return den


def _monic_associate(f, ring):
    qring = ring.clone(domain=ring.domain.get_field())
    f = _alpha_to_z(f.monic(), qring)
    f_ = f.clear_denoms()[1].set_ring(ring)

    return f_.primitive()[1]


def _leading_coeffs(f, U, gamma, lcfactors, A, D, denoms, divisors):
    ring = f.ring
    domain = ring.domain
    symbols = f.ring.symbols
    qring = ring.clone(symbols=(symbols[0], symbols[-1]), domain=domain.get_field())

    omega = domain(D * gamma.rep[0])

    denominators = [_denominator(u, qring) for u, _ in U]

    m = len(denoms)

    gcd = domain.gcd

    for i in xrange(m):
        pi = gcd(omega, divisors[i])
        divisors[i] //= pi
        if divisors[i] == 1:
            return None

    e = []

    for dj in denominators:
        ej = []

        for i in xrange(m):
            eji = 0
            g1 = gcd(dj, divisors[i])

            while g1 != 1:
                eji += 1
                dj = dj // g1
                g1 = gcd(dj, g1)

            ej.append(eji)

        e.append(ej)

    n = len(denominators)
    if any(sum([e[j][i] for j in xrange(n)]) != lcfactors[i][1] for i in xrange(m)):
        return None

    lcs = []
    for j in xrange(n):
        lj = ring.drop(0).mul([lcfactors[i][0]**e[j][i] for i in xrange(m)])
        lcs.append(lj)

    zring = qring.clone(domain=domain)
    U_ = [_alpha_to_z(u, qring).clear_denoms()[1].set_ring(zring) for u, _ in U]

    if omega == 1:
        for j in xrange(n):
            dj = lcs[j]
            djA = dj.evaluate(zip(dj.ring.gens[:-1], A)).drop(0)
            lcuj = U_[j].LC

            lcs[j] = dj.mul_ground(lcuj // djA)
    else:
        for j in xrange(n):
            dj = lcs[j]
            djA = dj.evaluate(zip(dj.ring.gens[:-1], A)).drop(0)
            lcuj = U_[j].LC
            d = gcd(djA, lcuj)

            lcs[j] = dj.mul_ground(lcuj // d)
            U_[j] = U_[j].mul_ground(djA // d)
            omega = (omega * d) // djA

        if omega == 1:
            return f_, lcs, U_
        else:
            lcs = [lc.mul_ground(omega) for lc in lcs]
            U_ = [u.mul_ground(omega) for u in U_]
            f_ = f.mul_ground(omega**(n-1))

    return f, lcs, U_


def _test_evaluation_points(f, gamma, lcfactors, D, A):
    ring = f.ring
    qring = ring.clone(domain=ring.domain.get_field())

    lc = ring.dmp_LC(f)
    if lc.evaluate(zip(lc.ring.gens, A)) == 0:
        return None

    fA = f.evaluate(zip(ring.gens[1:-1], A))
    if not fA.is_squarefree:
        return None

    omega = gamma * D
    denoms = []
    for l, _ in lcfactors:
        lA = l.evaluate(zip(l.ring.gens, A)) # in Q(alpha)
        denoms.append(_denominator(lA**(-1), qring))

    divisors = _distinct_prime_divisors(denoms, ring.domain)

    if divisors is None:
        return None
    elif any(omega % d == 0 for d in divisors):
        return None

    return fA, denoms, divisors


# TODO!
def _padic_lift(f, pfactors, l, p, B, lcs, minpoly):
    ring = f.ring
    LC = ring.dmp_LC
    x = ring.gens[0]

    h = [g + (l - LC(g))*x**g.degree() for g, l in zip(pfactors, lcs)]

    e = f - ring.mul(h).mul_ground(l) # mod minpoly

    P = p


def _div(f, g, minpoly, p):
    ring = f.ring

    rem = f
    deg = g.degree(0)
    lcinv, _, gcd = _gf_gcdex(ring.dmp_LC(g), minpoly, p)

#    minpoly mod p should be irreducible at this point
    if not gcd == 1:
        return None

    quotient = ring.zero

    while True:
        degrem = rem.degree(0)
        if degrem < deg:
            break
        quo = (lcinv * ring.dmp_LC(rem)).set_ring(ring).mul_monom((degrem - deg, 0))
        rem = _trunc(rem - g*quo, minpoly, p)
        quotient += quo

    return _trunc(quotient, minpoly, p), rem


def _extended_euclidean_algorithm(f, g, minpoly, p):
    ring = f.ring
    zero = ring.zero
    one = ring.one

    f = _trunc(f, minpoly, p)
    g = _trunc(g, minpoly, p)

    s0, s1 = zero, one
    t0, t1 = one, zero

    while g:
        quo, rem = _div(f, g, minpoly, p)
        f, g = g, rem
        s0, s1 = s1 - quo*s0, s0
        t0, t1 = t1 - quo*t0, t0

    lcfinv = _gf_gcdex(ring.dmp_LC(f), minpoly, p)[0].set_ring(ring)

    return ( _trunc(s1 * lcfinv, minpoly, p), _trunc(t1 * lcfinv, minpoly, p),
        _trunc(f * lcfinv, minpoly, p) )


def _diophantine_univariate(F, m, minpoly, p):
    if len(F) == 2:
        f, g = F
        s, t, _ = _extended_euclidean_algorithm(g, f, minpoly, p)

        s = s.mul_monom((m, 0))
        t = t.mul_monom((m, 0))

        q, s = _div(s, f, minpoly, p)

        t += q*g

        s = _trunc(s, minpoly, p)
        t = _trunc(t, minpoly, p)

        result = [s, t]
    else:
        ring = F[0].ring
        G = [F[-1]]

        for f in reversed(F[1:-1]):
            G.insert(0, f * G[0])

        S, T = [], [ring.one]

        for f, g in zip(F, G):
            t, s = _diophantine([g, f], T[-1], [], 0, minpoly, p)
            T.append(t)
            S.append(s)

        result, S = [], S + [T[-1]]

        for s, f in zip(S, F):
            r = _div(s.mul_monom((m, 0)), f, minpoly, p)[1]
            s = _trunc(r, minpoly, p)

            result.append(s)

    return result


def _diophantine(F, c, A, d, minpoly, p):
    ring = c.ring

    if not A:
        S = [ring.zero for _ in F]
        c = _trunc(c, minpoly, p)

        for (exp,), coeff in c.drop_to_ground(1).iterterms():
            T = _diophantine_univariate(F, exp, minpoly, p)

            for j, (s, t) in enumerate(zip(S, T)):
                S[j] = _trunc(s + t*coeff.set_ring(ring), minpoly, p)
    else:
        n = len(A)
        e = ring.mul(F)

        a, A = A[-1], A[:-1]
        B, G = [], []

        for f in F:
            B.append(e.quo(f))
            G.append(f.evaluate(n, a))

        C = c.evaluate(n, a)

        S = _diophantine(G, C, A, d, minpoly, p)
        S = [s.set_ring(ring) for s in S]

        for s, b in zip(S, B):
            c = c - s*b

        c = _trunc(c, minpoly, p)

        m = ring.gens[n] - a
        M = ring.one

        for k in xrange(d):
            if not c:
                break

            M = M * m
            C = ring.dmp_diff_eval_in(c, k + 1, a, n)

            if C:
                C = C.quo_ground(ring.domain.factorial(k + 1))
                T = _diophantine(G, C, A, d, minpoly, p)

                for i, t in enumerate(T):
                    T[i] = t.set_ring(ring) * M

                for i, (s, t) in enumerate(zip(S, T)):
                    S[i] = s + t

                for t, b in zip(T, B):
                    c = c - t * b

                c = _trunc(c, minpoly, p)

        S = [_trunc(s, minpoly, p) for s in S]

    return S


def _hensel_lift(f, H, LC, A, minpoly, p):
    ring = f.ring
    n = len(A)

    S = [f]
    H = list(H)

    for i, a in enumerate(reversed(A[1:])):
        s = S[0].evaluate(n - i, a)
        S.insert(0, _trunc(s, minpoly, p))

    d = max(f.degrees()[1:-1])

    for j, s, a in zip(xrange(1, n + 1), S, A):
        G = list(H)

        I, J = A[:j - 1], A[j:]

        Hring = f.ring
        for _ in xrange(j, n):
            Hring = Hring.drop(j + 1)

        x = Hring.gens[0]
        evalpoints = zip(LC[0].ring.gens[j:-1], J)

        for i, (h, lc) in enumerate(zip(H, LC)):
            if evalpoints:
                lc = lc.evaluate(evalpoints)
            lc = _trunc(lc, minpoly, p).set_ring(Hring)
            H[i] = h.set_ring(Hring) + (lc - h.LC)*x**h.degree()

        m = Hring.gens[j] - a
        M = Hring.one

        c = s - ring.mul(H)

        dj = s.degree(j)

        for k in xrange(dj):
            if not c:
                break

            M = M * m
            C = c.ring.dmp_diff_eval_in(c, k + 1, a, j)

            if C:
                C = C.quo_ground(ring.domain.factorial(k + 1)) # coeff of (x_{j-1} - a_{j-1})^(k + 1) in c
                T = _diophantine(G, C, I, d, minpoly, p)

                for i, (h, t) in enumerate(zip(H, T)):
                    H[i] = _trunc(h + t.set_ring(Hring)*M, minpoly, p)

                c = _trunc(s - ring.mul(H), minpoly, p)

    prod = ring.mul(H)
    if prod.rem(minpoly.set_ring(Hring)) != f:
        return None
    else:
        return H


def _sqf_p(f, minpoly, p):
    ring = f.ring
    lcinv, _, gcd = _gf_gcdex(ring.dmp_LC(f), minpoly, p)

    f = _trunc(f * lcinv.set_ring(ring), minpoly, p)

    if not f:
        return True
    else:
        return _euclidean_algorithm(f, _trunc(f.diff(0), minpoly, p), minpoly, p) == 1


def _test_prime(fA, minpoly, p, domain):
    if fA.LC % p == 0 or minpoly.LC % p == 0:
        return False
    if not _sqf_p(fA, minpoly, p):
        return False
    if not gf_irreducible_p(minpoly.to_dense(), p, domain):
        return False
    return True


# squarefree f with cont_x0(f) = 1
def _factor(f):
    ring = f.ring # Q(alpha)[x_0, ..., x_{n-1}]
    lcring = ring.drop(0)
    ground = ring.domain.domain
    n = ring.ngens

    if n == 1:
        lc, factors = f.factor_list()
        return (lc, [g for g, _ in factors])

    z = Dummy('z')
    qring = ring.clone(symbols=ring.symbols + (z,), domain=ground)
    lcqring = qring.drop(0)

    zring = qring.clone(domain=ground.get_ring())

    minpoly = _minpoly_from_dense(ring.domain.mod, zring.drop(*zring.gens[:-1]))
    f_ = _monic_associate(f, zring)
    D = minpoly.resultant(minpoly.diff(0))

#    # heuristic bound for p-adic lift
#    B = f_.max_norm()

    lc = zring.dmp_LC(f_)
    gamma, lcfactors = efactor(_z_to_alpha(lc, lcring)) # over QQ(alpha)[x_1, ..., x_n]

    # TODO: check if the computations of D_ and gamma_ are correct
    D_ = zring.domain.one
    gamma_ = gamma # in QQ(alpha)
    lcfactors_ = []

    for l, exp in lcfactors:
        den, l_ = _alpha_to_z(l, lcqring).clear_denoms() # l_ in QQ[x_1, ..., x_n, z], but coeffs in ZZ
        cont, l_ = l_.primitive()
        D_ *= den
        gamma_ *= cont
        lcfactors_.append((l_, exp)) # polyomials over QQ, allthough coeffs are in ZZ

    f_ = f_.mul_ground(D_)
    b = zring.dmp_zz_mignotte_bound(f_)
    p = nextprime(b)

    N = 0
    history = set([])
    tries = 5 # how big should this be?

    while True:
        for _ in xrange(tries):
            A = [random.randint(-N, N) for _ in xrange(n - 1)]

            if tuple(A) not in history:
                history.add(tuple(A))
            else:
                continue

            result = _test_evaluation_points(f_, gamma, lcfactors, D, A)
            if result is None:
                continue
            else:
                fA, denoms, divisors = result

            if any(denoms.count(denom) > 1 for denom in denoms):
                # TODO: check interval
                C = [random.randint(1, 3*(N + 1)) for _ in xrange(n - 1)]
                x = zring.gens[0]

                for i, ci in zip(xrange(1, n + 1), C):
                    xi = zring.gens[i]
                    f_ = f_.compose(xi, x + xi.mul_ground(ci))

                # TODO: check if this is still squarefree
                lc, factors = _factor(_z_to_alpha(f_, ring))

                for i, ci in zip(xrange(1, n + 1), C):
                    xi = zring.gens[i]
                    factors = [g.compose(xi, (xi - x).quo_ground(ci)) for g in factors]

                return (lc, factors)

            omega_, fAfactors = _z_to_alpha(fA, ring.drop(*ring.gens[1:])).factor_list() # factorization in Q(alpha)[x_0]
            if len(fAfactors) == 1:
                g = _z_to_alpha(f_, ring)
                return (f.LC, [g.monic()])

            result = _leading_coeffs(f_, fAfactors, gamma_, lcfactors_, A, D, denoms, divisors)
            if result is None:
                continue
            else:
                f_, lcs, fAfactors_ = result

            prod = ring.domain.domain.one
            for lc in lcs:
                prod *= lc.LC
            q = ground(prod, f_.LC)
            delta = ground.numer(q)
            l_ = ground.denom(q)

            f_ = f_.mul_ground(delta)

            while not _test_prime(fA, minpoly, p, zring.domain):
                p = nextprime(p)

            # what about l_ ?
            pfactors = _hensel_lift(f_, fAfactors_, lcs, A, minpoly, p)
            if pfactors is None:
                continue

#            factors = _padic_lift(f_, pfactors, l_, p, B, lcs, minpoly)
#            if factors is None:
#                B *= B
#                continue

            return (f.LC, [_z_to_alpha(g, ring).monic() for g in pfactors])

        N += 1


def _sort(factors, ring):
    densefactors = _sort_factors([(f.to_dense(), exp) for f, exp in factors])
    return [(ring.from_dense(f), exp) for f, exp in densefactors]


# output of the form (lc, [(poly1, exp1), ...])
def efactor(f):
    ring = f.ring

    assert ring.domain.is_Algebraic

    if f.is_ground:
        return (f[ring.zero_monom], [])

    n = ring.ngens

    if n == 1:
        return f.factor_list()
    else:
        cont, f = ring.dmp_primitive(f)
        if not cont.is_one:
            lccont, contfactors = efactor(cont)
            lc, factors = efactor(f)
            contfactors = [(g.set_ring(ring), exp) for g, exp in contfactors]
            return (lccont * lc, _sort(contfactors + factors, ring))

        # this is only correct because the content in x_0 is already divided out
        lc, sqflist = f.sqf_list()
        factors = []
        for g, exp in sqflist:
            lcg, gfactors = _factor(g)
            lc *= lcg
            factors = factors + [(gi, exp) for gi in gfactors]

        return (lc, _sort(factors, ring))
