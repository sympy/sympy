from sympy import (S, ceiling, floor, Dummy, Poly, resultant, quo, rem, Abs,
                   Add, cancel, monomials, apart_list, Lambda,
                   factorial, polygamma, RootSum, Mul, gruntz)

from sympy.core.sympify import sympify

from sympy.solvers.solvers import solve
from sympy.simplify.simplify import simplify


def integer_roots(p):
    r"""Compute integer roots of a polynomial.
    """
    intervals = p.intervals(eps=S(1)/4)
    roots = []

    for interval, mul in intervals:
        lo, hi = interval

        clo = ceiling(lo)
        fhi = floor(hi)

        if clo == fhi:
            roots.append((clo, mul))

    return roots


def dispersion(p, x):
    r"""Compute the 'dispersion' of a polynomial.

    For a polynomial `f(x)` with `deg f > 0` we define the dispersion
    as

    :math:`dis f(x) := max\{a \in Z^{+} \cup \{0\} | deg(gcd(f(x),f(x+a))) \geq 1\}`

    references: Abramov

    ..[1]: "On the Summation of Rational Functions"
    ..[2]: "Solutions of linear finite-differences equations with constant coefficients
            in the fields of rational functions"
    """
    a = Dummy("a")

    q = Poly(p.as_expr().subs(x, x+a), gens=[x])
    r = Poly(resultant(p, q), gens=[a])

    iroots = integer_roots(r)
    rmax = max([0] + [ r[0] for r in iroots ])

    return rmax


def step_1(p, q, l, x, dp, dq):
    # Split l**x * f1(x)/f2(x) into two parts h(x)*l**x and r(x)*l**x/f2(x)
    if dp >= dq:
        # p = h*q + r
        h = quo(p, q)
        r = rem(p, q)

        #print("Parts:")
        #print(h)
        #print(r)

        dr = r.degree()

        #print("Degree remainder:")
        #print(dr)

        if dr >= 0 and dr < dq:
            #print("SPLIT")

            s1 = fsum(h, l, x)
            s2 = fsum(r/q, l, x)

            return s1 + s2

    return None


def step_2(p, q, l, x, dp, dq):
    # Determine dispersion and alpha
    if dq == 0:
        alpha = 0
    else:
        ds = dispersion(q, x)
        alpha = ds

    #print("Dispersion (alpha):")
    #print(alpha)

    return alpha


def step_3(p, q, lam, x, dp, dq, alpha):
    # Find a degree bound
    if alpha == 0:
        if not Abs(lam) == 1:
            boundn = dp
        else:
            boundn = dp + 1

        g2 = 1
        p2 = 1

    # TODO: Really?
    #elif alpha > 1:
    elif alpha > 0:
        adds = []
        for i in xrange(0, alpha):
            tmp = p.shift(i)*lam**i / q.shift(i)
            adds.append(tmp)

        tmp = Add(*adds)
        tmp = cancel(tmp)

        s, t = tmp.as_numer_denom()
        s = Poly(s, gens=[x])
        t = Poly(t, gens=[x])

        tmp = cancel(t / t.shift(alpha))

        g2, unused = tmp.as_numer_denom()
        g2 = Poly(g2, gens=[x])

        p2 = cancel(p * g2 * g2.shift(1) / q)
        p2 = Poly(p2, gens=[x])
        #print("p2")
        #print(p2)

        l = p2.degree()
        m = g2.degree()

        #print("l:"+str(l))
        #print("m:"+str(m))

        if Abs(lam) == 1:
            if l > m-2:
                boundn = l - m + 1
            else:
                boundn = m
        else:
            boundn = l - m
    else:
        raise ValueError("What to do for alpha=1 ?")

    #print("g2:")
    #print(g2)

    #print("Degree bound:")
    #print(boundn)

    return boundn, g2, p2


def step_4(f, p, q, l, x, alpha, boundn, g2, p2):
    # Set up a candidate and try to solve for unknown coefficients
    from sympy.integrals.heurisch import _symbols
    ci = _symbols("c", boundn+1)
    mons = monomials([x], boundn)

    #print("********")
    #print(ci)
    #print(mons)

    g1 = Add(*[cii*mon for cii,mon in zip(ci,mons)])
    g1 = Poly(g1, gens=[x])

    #print("Candidate g1:")
    #print(g1)
    #print(g1.degree())

    if alpha == 0:
        # l*g1(x+1) - g1(x) = f(x)
        eqn = l * g1.shift(1).as_expr() - g1.as_expr() - f.as_expr()
    else:
        # p(x) = l*g1(x+1)*g2(x) - g2(x+1)*g1(x)
        eqn = l * g1.shift(1).as_expr()*g2.as_expr() - g2.shift(1).as_expr()*g1.as_expr() - p2.as_expr()
        #eqn = Poly(eqn, gens=[x]).as_expr()

    #print("Equation:")
    #pprint(eqn)

    sol = solve(eqn, ci, dict=True)

    # Hack for unhandy return format of solve:
    if type(sol) == list and len(sol) > 0:
        sol = sol[0]

    #print("Values:")
    #print(sol)

    return (g1, ci, sol)


def step_5(l, x, alpha, g1, g2, ci, sol):
    # Try to find a closed form in F
    # This is possible only if all ci got values
    # TODO: improve using linalg and solvability criteria
    unresolved = set(ci) - set(sol.keys())
    #print("Unresolved symbols:")
    #print(unresolved)

    if len(unresolved) <= 1:
        # System was solvable
        if alpha == 0:
            solution = l**x * g1.subs(sol)
        else:
            solution = l**x * g1.subs(sol) / g2.as_expr()
    else:
        if Abs(l) == 1:
            solution = None
        else:
            #solution = None
            ##print("Not summable in F or E ?")
            raise ValueError("Not summable in F or E   (1)")

    #print("Soluion:")
    #print(solution)

    return solution


def step_6(f, l, x):
    # Try for a closed form in an extension E of F
    cpfd = apart_list(f)

    #print("Common factor:")
    #print(cpfd[0])

    # Polynomial part
    f0 = cpfd[1]

    if f0.degree() >= 0:
        g0 = cpfd[0] * fsum(f0, l, x)
    else:
        g0 = 0

    #print("Polynomial part:")
    #print(g0)

    # Rational part
    # TODO: Handling of l wrong here!
    rp = []

    for r, nf, df, ex in cpfd[2]:
        an ,nu = nf.args
        alpha = an[0]
        term = nu * (-1)**(ex-1) / factorial(ex-1) * polygamma(ex-1, x-alpha)
        func = Lambda(alpha, term)
        rp.append(RootSum(r, func, auto=False))

    #print("Rational part:")
    #print(rp)

    return g0 + cpfd[0]*Add(*rp)


def preprocess(f, x):
    # Preprocess f to split off lambda
    f = sympify(f)

    #print(50*"=")
    #print("Entering preprocessor with:")
    #print("f: "+str(f))

    from sympy.integrals.heurisch import components
    monoms = components(f, x)

    #print("Monomials")
    #print(monoms)
    monoms = monoms - set([x])

    exppart = Mul(*list(monoms))
    exppart = simplify(exppart)

    #print("Exp-part")
    #print(exppart)

    apow = exppart.as_powers_dict()
    exps = apow.values()[0]

    tobase = cancel(exps / x)

    l = apow.keys()[0] ** tobase

    #print("lambda")
    #print(l)

    f = cancel(f / l**x)

    return f, l


def fsum(f, l, x):
    r"""Implements the algorithm `fsum` from the paper:
    'On Computing Closed Forms for Indefinite Summations'
    by Yiu-Kwong Man
    doi : "10.1006/jsco.1993.1053"
    """
    # The actual fsum algorithm
    #print(50*"=")
    #print("Entering fsum with:")
    #print("f: "+str(f))
    #print("l: "+str(l))
    #print(20*"-")

    # Split f
    #print(f)

    p, q = f.as_numer_denom()
    #print("p / q")
    #print(p)
    #print(q)

    p = Poly(p, gens=[x])
    q = Poly(q, gens=[x])

    dp = p.degree()
    dq = q.degree()

    #print("degrees")
    #print(dp)
    #print(dq)

    #
    result = step_1(p, q, l, x, dp, dq)

    if result is not None:
        return result

    #
    alpha = step_2(p, q, l, x, dp, dq)

    if dq > 0 and alpha == 0:
        if Abs(l) == 1:
            #print("Go to step S6")
            return step_6(f, l, x)
        else:
            ##print("Not summable in F and E ?")
            raise ValueError("Not summable in F and E   (2)")
            #return step_6(f, l, x)

    #
    bound, g2, p2 = step_3(p, q, l, x, dp, dq, alpha)
    #
    g1, coeffs, vals = step_4(f, p, q, l, x, alpha, bound, g2, p2)
    #
    solution = step_5(l, x, alpha, g1, g2, coeffs, vals)

    if solution is None:
        solution = step_6(f, l, x)

    # Do not call expand_func() here to resolve polygamma function if possible

    return solution


def ratsum(f, x):
    f, l = preprocess(f, x)
    return fsum(f, l, x)


def ratsum_def(f, bounds):
    x, a, b = bounds
    x = sympify(x)
    a = sympify(a)
    b = sympify(b)

    F = ratsum(f, x)

    # TODO: Better always use limits?
    if not a.is_finite:
        lo = gruntz(F, x, a)
    else:
        lo = F.subs(x, a)

    if not b.is_finite:
        hi = gruntz(F, x, b)
    else:
        hi = F.subs(x, b+1)

    Sn = hi - lo
    return Sn
