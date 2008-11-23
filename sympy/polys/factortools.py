"""General purpose factorization routines. """

from sympy.polys.integerpolys import zzx_degree, zzx_from_dict, zzx_factor, zzx_mul
from sympy.polys.polynomial import Poly, PolynomialError
from sympy.polys.monomial import monomial_div

from sympy.core import Integer, Rational, Symbol, sympify

def poly_factors(f, *symbols, **flags):
    """Factor polynomials over rationals.

       >>> from sympy import *
       >>> x, y = symbols("x y")

       >>> poly_factors(x**2 - y**2, x, y)
       (1, [(Poly(x - y, x, y), 1), (Poly(x + y, x, y), 1)])

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    denom, f = f.as_integer()

    if f.is_univariate:
        coeffs = map(int, f.iter_all_coeffs())
        content, factors = zzx_factor(coeffs)

        for i in xrange(len(factors)):
            factor, k = factors[i]
            n = zzx_degree(factor)

            terms = {}

            for j, coeff in enumerate(factor):
                if coeff != 0:
                    terms[(n-j,)] = Integer(coeff)

            factors[i] = Poly(terms, *f.symbols), k
    else:
        content, factors = kronecker_mv(f, **flags)

    return Rational(content, denom), factors

def factors(f, *symbols, **flags):
    """Factor polynomials over rationals.

       >>> from sympy import *
       >>> x, y = symbols("x y")

       >>> factors(x**2 - y**2, x, y)
       (1, [(x - y, 1), (x + y, 1)])

    """
    coeff, factors = poly_factors(f, *symbols, **flags)
    return coeff, [ (g.as_basic(), k) for g, k in factors ]

def factor(f, *symbols, **flags):
    """Factor polynomials over rationals.

       >>> from sympy import *
       >>> x, y = symbols("x y")

       >>> factor(x**2 - y**2)
       (x - y)*(x + y)

    """
    if not symbols and not isinstance(f, Poly):
        symbols = sympify(f).atoms(Symbol)

        if not symbols:
            return f

    coeff, factors = poly_factors(f, *symbols, **flags)

    result = 1 # XXX: don't include coeff in the leading factor

    for factor, k in factors:
        result *= factor.as_basic()**k

    return coeff * result

def kronecker_mv(f, **flags):
    """Kronecker method for Z[X] polynomials.

       NOTE: This function is very slow even on small input.
             Use debug=True flag to see its progress, if any.
    """
    symbols = f.symbols

    def mv_int_div(f, g):
        q = Poly((), *symbols)
        r = Poly((), *symbols)

        while not f.is_zero:
            lc_f, lc_g = f.LC, g.LC

            dv = lc_f % lc_g
            cf = lc_f / lc_g

            monom = monomial_div(f.LM, g.LM)

            if dv == 0 and monom is not None:
                q  = q.add_term(cf, monom)
                f -= g.mul_term(cf, monom)
            else:
                r = r.add_term(*f.LT)
                f = f.kill_lead_term()

        return q, r

    def combinations(lisp, m):
        def recursion(fa, lisp, m):
            if m == 0:
                 yield fa
            else:
                for i, fa2 in enumerate(lisp[0 : len(lisp) + 1 - m]):
                    for el in recursion(zzx_mul(fa2, fa), list(lisp[i + 1:]), m - 1):
                        yield el

        for i, fa in enumerate(lisp[0 : len(lisp) + 1 - m]):
            for el in recursion(fa, list(lisp[i + 1:]), m - 1):
                yield el

    debug = flags.get('debug', False)

    cont, f = f.as_primitive()
    N = len(symbols)

    max_exp = {}

    for v in symbols:
        max_exp[v] = 0

    for coeff, monom in f.iter_terms():
        for v, exp in zip(symbols, monom):
            if exp > max_exp[v]:
                max_exp[v] = exp

    symbols = sorted(symbols, reverse=True,
        key=lambda v: max_exp[v])

    f = Poly(f, *symbols)

    d = max_exp[symbols[0]] + 1

    terms, exps = {}, []

    for i in xrange(0, len(symbols)):
        exps.append(d**i)

    for coeff, monom in f.iter_terms():
        exp = 0

        for i, expi in enumerate(monom):
            exp += expi * exps[i]

        terms[exp] = int(coeff)

    g, factors = zzx_from_dict(terms), []

    try:
        for ff, k in zzx_factor(g)[1]:
            for i in xrange(0, k):
                factors.append(ff)
    except OverflowError:
        raise PolynomialError("input too large for multivariate Kronecker method")

    const, result, tested = 1, [], []

    if debug: print "KRONECKER-MV: Z[x] #factors = %i ..." % (len(factors))

    for k in range(1, len(factors)//2 + 1):
        for h in combinations(factors, k):
            if h in tested:
                continue

            n = zzx_degree(h)
            terms = {}

            for coeff in h:
                if not coeff:
                    n = n-1
                    continue
                else:
                    coeff = Integer(coeff)

                y_deg, n = n, n-1
                monom = [0] * N

                for i in xrange(N):
                    v_deg =  y_deg % d
                    y_deg = (y_deg - v_deg) // d
                    monom[i] = v_deg

                monom = tuple(monom)

                if terms.has_key(monom):
                    terms[monom] += coeff
                else:
                    terms[monom] = coeff

            cand = Poly(terms, *symbols)

            if cand.is_one:
                continue

            if cand.LC.is_negative:
                cand = -cand;

            q, r = mv_int_div(f, cand)

            if r.is_zero:
                if debug: print "KRONECKER-MV: Z[X] factor found %s" % cand
                result.append(cand)
                f = q
            else:
                tested.append(h)

            if f.is_constant:
                const, f = f.LC, Poly(1, *symbols)
                break

        if f.is_one:
            break

    if not f.is_one:
        if debug: print "KRONECKER-MV: Z[X] factor found %s" % f
        result.append(f)

    factors = {}

    for ff in result:
        if factors.has_key(ff):
            factors[ff] += 1
        else:
            factors[ff] = 1

    return cont*const, sorted(factors.items())

