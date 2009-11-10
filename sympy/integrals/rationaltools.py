"""This module implements tools for integrating rational functions. """

from sympy import S, Symbol, symbols, I, log, atan, Poly, \
    div, quo, resultant, roots, collect, solve, RootSum, Lambda

from sympy.polys.algorithms import poly_div, poly_gcd, \
    poly_gcdex, poly_sqf, poly_subresultants
from sympy.polys.rootfinding import number_of_real_roots

def ratint(f, x, **flags):
    """Performs indefinite integration of rational functions.

       Given a field K and a rational function f = p/q, where p and q
       are polynomials in K[x], returns a function g such that f = g'.

       >>> from sympy.integrals.rationaltools import ratint
       >>> from sympy.abc import x

       >>> ratint(36/(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2), x)
       -4*log(1 + x) + 4*log(-2 + x) - (6 + 12*x)/(1 - x**2)

       References
       ==========

       .. [Bro05] M. Bronstein, Symbolic Integration I: Transcendental
          Functions, Second Edition, Springer-Verlang, 2005, pp. 35-70

    """
    if type(f) is not tuple:
        p, q = f.as_numer_denom()
    else:
        p, q = f

    p, q = Poly(p, x), Poly(q, x)

    g = poly_gcd(p, q)

    p = poly_div(p, g)[0]
    q = poly_div(q, g)[0]

    result, p = poly_div(p, q)

    result = result.integrate(x).as_basic()

    if p.is_zero:
        return result

    g, h = ratint_ratpart(p, q, x)

    P, Q = h.as_numer_denom()
    q, r = poly_div(P, Q, x)

    result += g + q.integrate(x).as_basic()

    if not r.is_zero:
        symbol = flags.get('symbol', 't')

        if not isinstance(symbol, Symbol):
            t = Symbol(symbol, dummy=True)
        else:
            t = symbol

        L = ratint_logpart(r, Q, x, t)

        real = flags.get('real')

        if real is None:
            if type(f) is not tuple:
                atoms = f.atoms()
            else:
                p, q = f

                atoms = p.atoms() \
                      | q.atoms()

            for elt in atoms - set([x]):
                if not elt.is_real:
                    real = False
                    break
            else:
                real = True

        eps = S(0)

        if not real:
            for h, q in L:
                eps += RootSum(Lambda(t, t*log(h.as_basic())), q)
        else:
            for h, q in L:
                R = log_to_real(h, q, x, t)

                if R is not None:
                    eps += R
                else:
                    eps += RootSum(Lambda(t, t*log(h.as_basic())), q)

        result += eps

    return result

def ratint_ratpart(f, g, x):
    """Horowitz-Ostrogradsky algorithm.

       Given a field K and polynomials f and g in K[x], such that f and g
       are coprime and deg(f) < deg(g), returns fractions A and B in K(x),
       such that f/g = A' + B and B has square-free denominator.

    """
    f, g = Poly(f, x), Poly(g, x)

    u = poly_gcd(g, g.diff())
    v = poly_div(g, u)[0]

    n = u.degree - 1
    m = v.degree - 1
    d = g.degree

    A_coeff = [ Symbol('a' + str(n-i), dummy=True) for i in xrange(0, n+1) ]
    B_coeff = [ Symbol('b' + str(m-i), dummy=True) for i in xrange(0, m+1) ]

    symbols = A_coeff + B_coeff

    A = Poly(zip(A_coeff, xrange(n, -1, -1)), x)
    B = Poly(zip(B_coeff, xrange(m, -1, -1)), x)

    H = f - A.diff()*v + A*poly_div(u.diff()*v, u)[0] - B*u

    result = solve(H.coeffs, symbols)

    A = A.subs(result)
    B = B.subs(result)

    rat_part = Poly.cancel((A, u), x)
    log_part = Poly.cancel((B, v), x)

    return rat_part, log_part

def ratint_logpart(f, g, x, t=None):
    """Lazard-Rioboo-Trager algorithm.

       Given a field K and polynomials f and g in K[x], such that f and g
       are coprime, deg(f) < deg(g) and g is square-free, returns a list
       of tuples (s_i, q_i) of polynomials, for i = 1..n, such that s_i
       in K[t, x] and q_i in K[t], and:
                               ___    ___
                     d  f   d  \  `   \  `
                     -- - = --  )      )   a log(s_i(a, x))
                     dx g   dx /__,   /__,
                              i=1..n a | q_i(a) = 0

    """
    f, g = Poly(f, x), Poly(g, x)

    t = t or Symbol('t', dummy=True)
    a, b = g, f - g.diff().mul_term(t)

    res, R = poly_subresultants(a, b)
    Q = poly_sqf(Poly(res, t))

    R_map, H, i = {}, [], 1

    for r in R:
        R_map[r.degree] = r

    for q in Q:
        if q.degree > 0:
            _, q = q.as_primitive()

            if g.degree == i:
                H.append((g, q))
            else:
                h = R_map[i]
                A = poly_sqf(h.LC, t)

                for j in xrange(0, len(A)):
                    T = poly_gcd(A[j], q)**(j+1)
                    h = poly_div(h, Poly(T, x))[0]

                # NOTE: h.LC is always invertible in K[t]
                inv, coeffs = Poly(h.LC, t).invert(q), [S(1)]

                for coeff in h.coeffs[1:]:
                    T = poly_div(inv*coeff, q)[1]
                    coeffs.append(T.as_basic())

                h = Poly(zip(coeffs, h.monoms), x)

                H.append((h, q))

        i += 1

    return H

def log_to_atan(f, g):
    """Convert complex logarithms to real arctangents.

       Given a real field K and polynomials f and g in K[x], with g != 0,
       returns a sum h of arctangents of polynomials in K[x], such that:

                       df   d         f + I g
                       -- = -- I log( ------- )
                       dx   dx        f - I g

    """
    if f.degree < g.degree:
        f, g = -g, f

    p, q = poly_div(f, g)

    if q.is_zero:
        return 2*atan(p.as_basic())
    else:
        s, t, h = poly_gcdex(g, -f)
        A = 2*atan(quo(f*s+g*t, h))

        return A + log_to_atan(s, t)

def log_to_real(h, q, x, t):
    """Convert complex logarithms to real functions.

       Given real field K and polynomials h in K[t,x] and q in K[t],
       returns real function f such that:
                              ___
                      df   d  \  `
                      -- = --  )  a log(h(a, x))
                      dx   dx /__,
                             a | q(a) = 0

    """
    u, v = symbols('u v')

    H = h.subs({t:u+I*v}).as_basic().expand()
    Q = q.subs({t:u+I*v}).as_basic().expand()

    H_map = collect(H, I, evaluate=False)
    Q_map = collect(Q, I, evaluate=False)

    a, b = H_map.get(S(1), S(0)), H_map.get(I, S(0))
    c, d = Q_map.get(S(1), S(0)), Q_map.get(I, S(0))

    R = Poly(resultant(c, d, v), u)

    R_u = roots(R, domain='R')

    if len(R_u) != number_of_real_roots(R):
        return None

    result = S(0)

    for r_u in R_u.iterkeys():
        C = Poly(c.subs({u:r_u}), v)
        R_v = roots(C, domain='R')

        if len(R_v) != number_of_real_roots(C):
            return None

        for r_v in R_v:
            if not r_v.is_positive:
                continue

            D = d.subs({u:r_u, v:r_v})

            if D.evalf(chop=True) != 0:
                continue

            A = Poly(a.subs({u:r_u, v:r_v}), x)
            B = Poly(b.subs({u:r_u, v:r_v}), x)

            AB = (A**2 + B**2).as_basic()

            result += r_u*log(AB) + r_v*log_to_atan(A, B)

    R_q = roots(q, domain='R')

    if len(R_q) != number_of_real_roots(q):
        return None

    for r in R_q.iterkeys():
        result += r*log(h.subs(t, r).as_basic())

    return result

