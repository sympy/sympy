"""Algorithms for partial fraction decomposition of rational functions. """

from sympy.core import S
from sympy.polys import Poly
from sympy.solvers import solve
from sympy.utilities import numbered_symbols, take

def apart_undetermined_coeffs(P, Q):
    """Partial fractions using method of undetermined coefficients. """
    common, P, Q = P.cancel(Q)
    poly_part, P = P.div(Q)

    _, factors = Q.factor_list()

    X = numbered_symbols(dummy=True)
    partial, symbols = [], []

    for f, k in factors:
        n, q = f.degree(), Q

        for i in xrange(1, k+1):
            coeffs, q = take(X, n), q.exquo(f)
            partial.append((coeffs, q, f, i))
            symbols.extend(coeffs)

    dom = Q.get_domain().inject(*symbols)
    F = Poly(0, Q.gen, domain=dom)

    for i, (coeffs, q, f, k) in enumerate(partial):
        h = Poly(coeffs, Q.gen, domain=dom)
        partial[i] = (h, f, k)
        q = q.set_domain(dom)
        F += h*q

    system, result = [], S(0)

    for (k,), coeff in F.terms():
        system.append(coeff - P.nth(k))

    solution = solve(system, symbols)

    for h, f, k in partial:
        h = h.as_basic().subs(solution)
        result += h/f.as_basic()**k

    return common*(poly_part.as_basic() + result)
