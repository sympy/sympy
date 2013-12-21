from __future__ import print_function, division

from sympy.core import Add, Mul
from sympy.core.sympify import sympify
from sympy.core.function import Lambda
from sympy.utilities import flatten
from sympy.simplify.simplify import simplify
from sympy.functions.elementary.complexes import Abs
from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.special.gamma_functions import polygamma
from sympy.polys.monomials import itermonomials
from sympy.polys.polytools import Poly, quo, rem, cancel
from sympy.polys.partfrac import apart_list
from sympy.polys.rootoftools import RootSum
from sympy.concrete.dispersion import dispersion
from sympy.matrices import zeros
from sympy.solvers.linear import solve_general_linear
from sympy.series import gruntz


def split_division(p, q, l, x):
    r"""
    """
    # Split l**x * f1(x)/f2(x) into two parts h(x)*l**x and r(x)*l**x/f2(x)
    dp = p.degree()
    dq = q.degree()

    if dp >= dq:
        # p = h*q + r
        h = quo(p, q)
        r = rem(p, q)
        dr = r.degree()

        if dr >= 0 and dr < dq:
            s1 = fsum(h, l, x)
            s2 = fsum(r/q, l, x)
            return s1 + s2

    return None


def compute_dispersion(p, q, l, x):
    r"""
    """
    # Determine dispersion and alpha
    dq = q.degree()

    if dq == 0:
        alpha = 0
    else:
        ds = dispersion(q)
        alpha = ds
    return alpha


def degree_bound(p, q, lam, x, alpha):
    r"""
    """
    # Find a degree bound
    dp = p.degree()

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
        for i in range(0, alpha):
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

        l = p2.degree()
        m = g2.degree()

        if Abs(lam) == 1:
            if l > m-2:
                boundn = l - m + 1
            else:
                boundn = m
        else:
            boundn = l - m
    else:
        raise ValueError("What to do for alpha=1 ?")

    return boundn, g2, p2


def find_relation(f, p, q, l, x, alpha, boundn, g2, p2):
    r"""
    """
    # Set up a candidate and try to solve for unknown coefficients
    from sympy.integrals.heurisch import _symbols
    ci = _symbols("c", boundn+1)
    mons = itermonomials([x], boundn)

    g1 = Add(*[cii*mon for cii,mon in zip(ci,mons)])
    g1 = Poly(g1, gens=[x])

    if alpha == 0:
        # l*g1(x+1) - g1(x) = f(x)
        eqn = l * g1.shift(1).as_expr() - g1.as_expr() - f.as_expr()
    else:
        # p(x) = l*g1(x+1)*g2(x) - g2(x+1)*g1(x)
        eqn = l * g1.shift(1).as_expr()*g2.as_expr() - g2.shift(1).as_expr()*g1.as_expr() - p2.as_expr()

    P = Poly(eqn, gens=[x])

    M = zeros(len(mons), len(ci))
    rhs = zeros(len(mons), 1)

    for r, mon in enumerate(mons):
        cmon = P.coeff_monomial(mon)
        cmonp = Poly(cmon, gens=ci)
        # Fill in matrix elements
        for c, coe in enumerate(ci):
            xi = cmonp.coeff_monomial(coe)
            M[r,c] = xi
        # Put constants in the RHS
        chi = cmonp.coeff_monomial(1)
        rhs[r,0] = -chi

    try:
        solution = solve_general_linear(M, rhs)
    except ValueError:
        solution = (None, None)

    return (g1, ci, solution)


def closed_form_base_field(l, x, alpha, g1, g2, ci, linsol):
    # Try to find a closed form in F
    # This is possible only if all ci got values
    sol, params = linsol
    sol, params = flatten(sol.tolist()), flatten(params.tolist())
    vals = [(var,val) for var, val in zip(ci, sol)]

    if params is not None and len(params) <= 1:
        # System was solvable
        if alpha == 0:
            solution = l**x * g1.subs(vals)
        else:
            solution = l**x * g1.subs(vals) / g2.as_expr()
        # Remove constant of summation
        if len(params) == 1:
            solution = solution.subs(params[0], 0)
    else:
        if Abs(l) == 1:
            solution = None
        else:
            raise ValueError("Not summable in F or E   (1)")

    return solution


def closed_form_field_extension(f, l, x):
    # Try for a closed form in an extension E of F
    cpfd = apart_list(f)

    # Polynomial part
    f0 = cpfd[1]

    if f0.degree() >= 0:
        g0 = cpfd[0] * fsum(f0, l, x)
    else:
        g0 = 0

    # Rational part
    # TODO: Handling of l wrong here!
    rp = []

    for r, nf, df, ex in cpfd[2]:
        an ,nu = nf.args
        alpha = an[0]
        term = nu * (-1)**(ex-1) / factorial(ex-1) * polygamma(ex-1, x-alpha)
        func = Lambda(alpha, term)
        rp.append(RootSum(r, func, auto=False))

    return g0 + cpfd[0]*Add(*rp)


def preprocess(f, x):
    # Preprocess f to split off lambda
    f = sympify(f)

    from sympy.integrals.heurisch import components
    monoms = components(f, x)
    monoms = monoms - set([x])

    exppart = Mul(*list(monoms))
    exppart = simplify(exppart)

    apow = exppart.as_powers_dict()
    exps = list(apow.values())[0]

    tobase = cancel(exps / x)

    l = list(apow.keys())[0] ** tobase

    if x in l.free_symbols:
        raise ValueError("Base of exponential part not free of "+str(x)+".")

    f = cancel(f / l**x)

    return f, l


def fsum(f, l, x):
    r"""Implements the algorithm `fsum`.

    References
    ==========

    .. [1] Yiu-Kwong Man, "On Computing Closed Forms for Indefinite Summations",
           Journal of Symbolic Computation, Volume 16, Issue 4, 1993, Pages 355-376
           http://www.sciencedirect.com/science/article/pii/S0747717183710539
    .. [2] Yiu-Kwong Man and Francis J. Wright, "Fast Polynomial Dispersion
           Computation and its Application to Indefinite Summation",
           Proceedings of the International Symposium on Symbolic and
           Algebraic Computation, 1994, Pages 175-180
           http://doi.acm.org/10.1145/190347.190413
    """
    # The actual fsum algorithm
    p, q = f.as_numer_denom()
    p = Poly(p, gens=[x])
    q = Poly(q, gens=[x])

    result = split_division(p, q, l, x)

    if result is not None:
        return result

    dq = q.degree()
    alpha = compute_dispersion(p, q, l, x)

    if dq > 0 and alpha == 0:
        if Abs(l) == 1:
            return closed_form_field_extension(f, l, x)
        else:
            #print("Not summable in F and E ?")
            raise ValueError("Not summable in F and E   (2)")
            #return closed_form_field_extension(f, l, x)

    bound, g2, p2 = degree_bound(p, q, l, x, alpha)
    g1, coeffs, vals = find_relation(f, p, q, l, x, alpha, bound, g2, p2)
    solution = closed_form_base_field(l, x, alpha, g1, g2, coeffs, vals)

    if solution is None:
        solution = closed_form_field_extension(f, l, x)

    # Do not call expand_func() here to resolve polygamma function if possible

    return solution


def ratsum(f, x):
    r"""

    References
    ==========

    .. [1] Yiu-Kwong Man, "On Computing Closed Forms for Indefinite Summations",
           Journal of Symbolic Computation, Volume 16, Issue 4, 1993, Pages 355-376
           http://www.sciencedirect.com/science/article/pii/S0747717183710539
    .. [2] Yiu-Kwong Man and Francis J. Wright, "Fast Polynomial Dispersion
           Computation and its Application to Indefinite Summation",
           Proceedings of the International Symposium on Symbolic and
           Algebraic Computation, 1994, Pages 175-180
           http://doi.acm.org/10.1145/190347.190413
    """
    f, l = preprocess(f, x)
    return fsum(f, l, x)


def ratsum_def(f, bounds):
    r"""

    References
    ==========

    .. [1] Yiu-Kwong Man, "On Computing Closed Forms for Indefinite Summations",
           Journal of Symbolic Computation, Volume 16, Issue 4, 1993, Pages 355-376
           http://www.sciencedirect.com/science/article/pii/S0747717183710539
    .. [2] Yiu-Kwong Man and Francis J. Wright, "Fast Polynomial Dispersion
           Computation and its Application to Indefinite Summation",
           Proceedings of the International Symposium on Symbolic and
           Algebraic Computation, 1994, Pages 175-180
           http://doi.acm.org/10.1145/190347.190413
    """
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
