"""Polynomial division algorithms for use with class Polynomial"""

from sympy.polynomials.base import *

def div(f, g, var=None, order=None, coeff=None):
    """Devides f by the polynomials in g, returns quotients and remainder.
    """
    
    if not isinstance(g, list):
        g = [g]
    # Only f is checked, the rest is assumed to match.
    if not isinstance(f, Polynomial):
        f = sympify(f)
        g = map(lambda x: sympify(x), g)
        if isinstance(var, Symbol):
            var = [var]
        if var is None:
            var = merge_var(f.atoms(type=Symbol),
                            *[g_i.atoms(type=Symbol) for g_i in g])
        f = Polynomial(f, var=var, order=order)
        g = map(lambda x: Polynomial(x, var=var, order=order), g)
    
    # Begin computation.
    r = Polynomial(S.Zero, var=f.var, order=f.order)
    q = []
    for i in range(0,len(g)):
        q.append(r)

    while f.sympy_expr is not S.Zero:    # f != 0
        for g_i in g:
            if g_i.sympy_expr is S.Zero: # Avoid division by 0.
                continue
            # Check if leading term of f is divisible by that of g_i.
            td = term_div(f.coeffs[0], g_i.coeffs[0])
            if (coeff != 'int' or isinstance(td[0], Integer)) \
               and all([e.is_nonnegative for e in td[1:]]):
                quot = Polynomial(coeffs=(td,), var=f.var, order=f.order)
                q[g.index(g_i)] += quot
                f -= quot*g_i
                break
        else: # No division occured, add the leading term to remainder.
            lt = f.leading_term()
            r += lt
            f -= lt

    if len(q) == 1:
        return q[0], r
    else:
        return q, r


def gcd(f, g, var=None, order=None, coeff=None):
    """Greatest common divisor.
    """

    # Check if f is a Polynomial already, g is assumed to match.
    if not isinstance(f, Polynomial):
        f = sympify(f)
        g = sympify(g)
        if var is None:
            var = merge_var(f.atoms(type=Symbol), g.atoms(type=Symbol))
        f = Polynomial(f, var=var, order=order)
        g = Polynomial(g, var=var, order=order)

    # Check if we need to keep an integer factor.
    if coeff == 'int':
        # TODO: Check for coefficient type?
        cf, f = f.as_primitive()
        cg, g = g.as_primitive()
        c = Integer(numbers.gcd(int(cf), int(cg)))
    else:
        c = S.One

    if len(f.var) == 0: # Constant result.
        return Polynomial(c, var=f.var, order=f.order)
    elif len(f.var) == 1: # Use euclidean algorithm.
        while g.sympy_expr is not S.Zero:
            lc, g = g.as_monic()
            f, g = g, div(f, g)[-1]
    else: # Use lcm and product to get multivariate gcd.
        l = lcm(f, g)
        q, r = div(f*g, l)
        assert r.sympy_expr is S.Zero
        lc, f = q.as_monic()

    return Polynomial(coeffs=tuple([(c*t[0],) + t[1:] for t in f.coeffs]),
                      var=f.var, order=f.order)


def lcm(f, g, var=None, order=None, coeff=None):
    """Least common multiple"""

    # Check if f is a Polynomial already, g is assumed to match.
    if not isinstance(f, Polynomial):
        f = sympify(f)
        g = sympify(g)
        if var is None:
            var = merge_var(f.atoms(type=Symbol), g.atoms(type=Symbol))
        f = Polynomial(f, var=var, order=order)
        g = Polynomial(g, var=var, order=order)

    # Check if we need to keep an integer factor.
    if coeff == 'int':
        # TODO: Check for coefficient type?
        cf, f = f.as_primitive()
        cg, g = g.as_primitive()
        cf, cg = int(cf), int(cg)
        c = Integer(cf*cg/numbers.gcd(cf, cg))
    else:
        c = S.One

    if len(f.var) == 0: # Constant result.
        return Polynomial(c, var=f.var, order=f.order)
    elif len(f.var) == 1: # Use gcd to get univariate lcm.
        gcd_temp = gcd(f, g)
        q, r = div(f*g, gcd_temp)
        assert r.sympy_expr is S.Zero
        lc, f = q.as_monic()
    else:
        # Compute a lexicographic Groebner base of the ideal generated
        # by t*f and (t-1)*g, with unrelated t.
        from sympy.polynomials import groebner_

        t = Symbol('t', dummy=True)
        var = [t] + f.var
        G = groebner_.groebner([Polynomial(t*f.sympy_expr,
                                           var=var, order='1-el'),
                                Polynomial((t-1)*g.sympy_expr,
                                           var=var, order='1-el')],
                               reduced=True)
        # Now intersect this result with the polynomial ring in the
        # var in `var', that is, eliminate t.
        I = filter(lambda p: t not in p.sympy_expr.atoms(type=Symbol), G)
        # The intersection should be a principal ideal, that is generated
        # by a single polynomial.
        if not len(I) == 1:
            raise PolynomialException("No single generator for intersection.")
        f = Polynomial(I[0].sympy_expr, var=f.var, order=f.order)

    return Polynomial(coeffs=tuple([(c*t[0],) + t[1:] for t in f.coeffs]),
                      var=f.var, order=f.order)

    
