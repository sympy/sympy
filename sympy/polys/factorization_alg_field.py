from __future__ import print_function

from sympy import Dummy
import random

from sympy.ntheory import nextprime
from sympy.integrals.heurisch import _symbols
from sympy.polys.galoistools import gf_irreducible_p
from sympy.polys.modulargcd import _trunc, _gf_gcdex, _minpoly_from_dense, _euclidean_algorithm
from sympy.polys.polyclasses import ANP
from sympy.polys.polyerrors import UnluckyLeadingCoefficient
from sympy.polys.polyutils import _sort_factors
from sympy.polys.rings import PolyRing
from sympy.polys.solvers import solve_lin_sys


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
                    f_[m] = coeff[i]

    return f_


def _z_to_alpha(f, ring):
    domain = ring.domain

    f_ = ring.zero
    for monom, coeff in f.iterterms():
        m = monom[:-1]
        c = domain([domain.domain(coeff)] + [0]*monom[-1])

        if m not in f_:
            f_[m] = c
        else:
            f_[m] += c

    return f_


def _distinct_prime_divisors(S, domain):
    gcd = domain.gcd
    divisors = []

    for i, s in enumerate(S):
        divisors.append(s)

        for j in xrange(i):
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


def _denominator(f):

    ring = f.ring.domain.get_ring()
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
    gcd = domain.gcd

    U = [_alpha_to_z(u, qring) for u, _ in U]
    denominators = [_denominator(u) for u in U]

    omega = D * gamma

    m = len(denoms)

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

    lcring = ring.drop(0)
    lcs = []
    for j in xrange(n):
        lj = lcring.mul([lcfactors[i][0]**e[j][i] for i in xrange(m)])
        lcs.append(lj)

    zring = qring.clone(domain=domain)

    for j in xrange(n):
        lj = lcs[j]
        dj = denominators[j]
        ljA = lj.evaluate(zip(lcring.gens, A))

        lcs[j] = lj.mul_ground(dj)
        U[j] = U[j].mul_ground(dj).set_ring(zring) * ljA.set_ring(zring)

        if omega == 1:
            f = f.mul_ground(dj)
        else:
            d = gcd(omega, dj)
            f = f.mul_ground(dj // d)

    if omega != 1:
        lcs[0] = lcs[0].mul_ground(omega)
        U[0] = U[0].mul_ground(omega)

    return f, lcs, U


def _test_evaluation_points(f, gamma, lcfactors, A, D):
    ring = f.ring
    qring = ring.clone(domain=ring.domain.get_field())

    fA = f.evaluate(zip(ring.gens[1:-1], A))

    if fA.degree() < f.degree():
        return None

    if not fA.is_squarefree:
        return None

    omega = gamma * D
    denoms = []
    for l, _ in lcfactors:
        lA = l.evaluate(zip(l.ring.gens, A)) # in Q(alpha)
        denoms.append(_denominator(_alpha_to_z(lA**(-1), qring)))

    if any(denoms.count(denom) > 1 for denom in denoms):
        raise UnluckyLeadingCoefficient

    divisors = _distinct_prime_divisors(denoms, ring.domain)

    if divisors is None:
        return None
    elif any(omega % d == 0 for d in divisors):
        return None

    return fA, denoms, divisors


def _subs_ground(f, A):
    f_ = f.ring.zero

    for monom, coeff in f.iterterms():
        if coeff.subs(A):
            f_[monom] = coeff.subs(A)

    return f_


def _choose_particular_solution(solution, ring):
    domain = ring.domain
    gens = list(ring.gens)
    sol = {}

    for k, v in solution.items():
        sol[k] = v.coeff(1)
        gens.remove(k)

    for k in gens:
        sol[k] = domain.zero

    return sol


def _padic_lift(f, pfactors, lcs, B, minpoly, p):
    ring = f.ring
    domain = ring.domain
    LC = ring.dmp_LC

    coeffs = []
    for i, g in enumerate(pfactors):
        coeffs += _symbols('c%i' % i, len(g))

    coeffring = PolyRing(coeffs, domain)
    ring_ = ring.clone(domain=coeffring)

    S = []
    k = 0
    for g in pfactors:
        s = ring_.zero
        t = len(g)
        for i, monom in zip(xrange(k, k+t), g.itermonoms()):
            s[monom] = coeffring.gens[i]
        S.append(s)
        k += t

    m = minpoly.set_ring(ring_)
    f = f.set_ring(ring_)
    x = ring_.gens[0]
    H = [g.set_ring(ring_) + (li - LC(g)).set_ring(ring_)*x**g.degree() for
            g, li in zip(pfactors, lcs)]

    prod = ring_.mul(H)
    e = (f - prod).rem(m)

    P = p
    while e and P < 2*B:
        poly = e // P

        for s, h in zip(S, H):
            poly -= prod.quo(h)*s

        poly = _trunc(poly, m, P)

        iszerofunc = lambda x: not (x % p)
        scalefunc = lambda x, _, scale: (x * domain.invert(scale, P)) % P
        elimfunc = lambda x, y, scale: (x - scale*y) % P

        solution = solve_lin_sys(poly.coeffs(), coeffring, iszerofunc=iszerofunc,
            scalefunc=scalefunc, elimfunc=elimfunc)

        if solution is None:
            return None
        else:
            for k, v in solution.items():
                solution[k] = v.trunc_ground(P)

        solution = _choose_particular_solution(solution, coeffring)
        subs = solution.items()

        H = [h + _subs_ground(s, subs).mul_ground(P) for h, s in zip(H, S)]
        P = P**2
        prod = ring_.mul(H)
        e = (f - prod).rem(m)

    if e == 0:
        return [h.set_ring(ring) for h in H]
    else:
        return None


def _div(f, g, minpoly, p):
    ring = f.ring

    rem = f
    deg = g.degree(0)
    lcinv, _, gcd = _gf_gcdex(ring.dmp_LC(g), minpoly, p)

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
        result = _div(f, g, minpoly, p)
        if result is None:
            return None
        else:
            quo, rem = result
        f, g = g, rem
        s0, s1 = s1 - quo*s0, s0
        t0, t1 = t1 - quo*t0, t0

    lcfinv = _gf_gcdex(ring.dmp_LC(f), minpoly, p)[0].set_ring(ring)

    return ( _trunc(s1 * lcfinv, minpoly, p), _trunc(t1 * lcfinv, minpoly, p),
        _trunc(f * lcfinv, minpoly, p) )


def _diophantine_univariate(F, m, minpoly, p):
    if len(F) == 2:
        f, g = F
        result = _extended_euclidean_algorithm(g, f, minpoly, p)
        if result is None:
            return None
        else:
            s, t, _ = result

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
            result = _diophantine([g, f], T[-1], [], 0, minpoly, p)
            if result is None:
                return None
            else:
                t, s = result
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
            if T is None:
                return None

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
        if S is None:
            return None
        else:
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
                if T is None:
                    return None

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
            H[i] = h.set_ring(Hring) + (lc - h.ring.dmp_LC(h).set_ring(Hring))*x**h.degree()

        m = Hring.gens[j] - a
        M = Hring.one

        c = _trunc(s - ring.mul(H), minpoly, p)

        dj = s.degree(j)

        for k in xrange(dj):
            if not c:
                break

            M = M * m
            C = c.ring.dmp_diff_eval_in(c, k + 1, a, j)

            if C:
                C = C.quo_ground(ring.domain.factorial(k + 1)) # coeff of (x_{j-1} - a_{j-1})^(k + 1) in c
                T = _diophantine(G, C, I, d, minpoly, p)
                if T is None:
                    return None

                for i, (h, t) in enumerate(zip(H, T)):
                    H[i] = _trunc(h + t.set_ring(Hring)*M, minpoly, p)

                c = _trunc(s - ring.mul(H), minpoly, p)

    prod = ring.mul(H)
    if _trunc(prod, minpoly, p) != f.trunc_ground(p):
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


def _test_prime(fA, D, minpoly, p, domain):
    if fA.LC % p == 0 or minpoly.LC % p == 0:
        return False
    if not _sqf_p(fA, minpoly, p):
        return False
    if D % p == 0:
        return False

    return True


# squarefree f with cont_x0(f) = 1
def _factor(f, save):
    ring = f.ring # Q(alpha)[x_0, ..., x_{n-1}]
    lcring = ring.drop(0)
    uniring = ring.drop(*ring.gens[1:])
    ground = ring.domain.domain
    n = ring.ngens

    z = Dummy('z')
    qring = ring.clone(symbols=ring.symbols + (z,), domain=ground)
    lcqring = qring.drop(0)

    groundring = ground.get_ring()
    zring = qring.clone(domain=groundring)
    lczring = zring.drop(0)

    minpoly = _minpoly_from_dense(ring.domain.mod, zring.drop(*zring.gens[:-1]))
    f_ = _monic_associate(f, zring)

    if save is True:
        D = minpoly.resultant(minpoly.diff(0))
    else:
        D = groundring.one

    # heuristic bound for p-adic lift
    B = (f_.max_norm() + 1)*D

    lc = zring.dmp_LC(f_)
    gamma, lcfactors = efactor(_z_to_alpha(lc, lcring)) # over QQ(alpha)[x_1, ..., x_n]

    gamma = ground.convert(gamma)
    D_ = ground.denom(gamma)
    gamma_ = ground.numer(gamma)
    lcfactors_ = []

    for l, exp in lcfactors:
        den, l_ = _alpha_to_z(l, lcqring).clear_denoms() # l_ in QQ[x_1, ..., x_n, z], but coeffs in ZZ
        cont, l_ = l_.set_ring(lczring).primitive()
        D_ *= den**exp
        gamma_ *= cont**exp
        lcfactors_.append((l_, exp))

    f_ = f_.mul_ground(D_)
    p = 2

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

            try:
                result = _test_evaluation_points(f_, gamma_, lcfactors, A, D)
            except UnluckyLeadingCoefficient:
                # TODO: check interval
                C = [random.randint(1, 3*(N + 1)) for _ in xrange(n - 1)]
                gens = zring.gens
                x = gens[0]

                for i, ci in zip(xrange(1, n + 1), C):
                    xi = gens[i]
                    f_ = f_.compose(xi, x + xi.mul_ground(ci))

                lc, factors = _factor(_z_to_alpha(f_, ring), save)
                gens = factors[0].ring.gens
                x = gens[0]

                for i, ci in zip(xrange(1, n + 1), C):
                    xi = gens[i]
                    factors = [g.compose(xi, (xi - x).quo_ground(ci)) for g in factors]

                return (lc, factors)

            if result is None:
                continue
            else:
                fA, denoms, divisors = result

            _, fAfactors = uniring.dup_ext_factor(_z_to_alpha(fA, uniring))
            if len(fAfactors) == 1:
                g = _z_to_alpha(f_, ring)
                return (f.LC, [g.monic()])

            result = _leading_coeffs(f_, fAfactors, gamma_, lcfactors_, A, D, denoms, divisors)
            if result is None:
                continue
            else:
                f_, lcs, fAfactors_ = result

            prod = groundring.one
            for lc in lcs:
                prod *= lc.LC
            delta = ground.numer(ground(prod, f_.LC))

            f_ = f_.mul_ground(delta)

            while not _test_prime(fA, D, minpoly, p, zring.domain):
                p = nextprime(p)

            pfactors = _hensel_lift(f_, fAfactors_, lcs, A, minpoly, p)
            if pfactors is None:
                p = nextprime(p)
                f_ = f_.primitive()[1]
                continue

            factors = _padic_lift(f_, pfactors, lcs, B, minpoly, p)
            if factors is None:
                p = nextprime(p)
                f_ = f_.primitive()[1]
                B *= B
                continue

            return (f.LC, [_z_to_alpha(g.primitive()[1], ring).monic() for g in factors])

        N += 1


def _sort(factors, ring):
    densefactors = _sort_factors([(f.to_dense(), exp) for f, exp in factors])
    return [(ring.from_dense(f), exp) for f, exp in densefactors]


# output of the form (lc, [(poly1, exp1), ...])
def efactor(f, save=True):
    ring = f.ring

    assert ring.domain.is_Algebraic

    if f.is_ground:
        return (f[ring.zero_monom], [])

    n = ring.ngens

    if n == 1:
        return ring.dup_ext_factor(f)
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
            lcg, gfactors = _factor(g, save)
            lc *= lcg
            factors = factors + [(gi, exp) for gi in gfactors]

        return (lc, _sort(factors, ring))
