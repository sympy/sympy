"""Various algorithms for the factorization of polynomials."""

import random

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import div_, gcd_, roots_, sqf_

def uv_int(f):
    """Find the factorization of an univariate integer polynomial.

    Using square-free factorization and Kronecker's algorithm."""
    a = sqf_.uv_int(f)
    result = []
    for i, p in enumerate(a):
        if p.coeffs[0][1] is not S.Zero: # Filter out constant factors
            # Check for rational roots first:
            rr = roots_.rat_roots(p)
            # In a square-free polynomial, the roots appear only once.
            for root in rr:
                pp = Polynomial(coeffs=((Rational(root.q), S.One),
                                       (Rational(-root.p), S.Zero)),
                                var=f.var, order=f.order)
                result += [pp]*(i+1)
                # TODO: remove assertion, speed up!
                q, r = div_.mv_int(p, pp)
                assert r.sympy_expr is S.Zero
                p = q[0] # q is a list!
            if p.coeffs[0][1] is not S.Zero: # Filter out constant factors
                # Then try the rest with the kronecker algorithm:
                for pp in kronecker(p):
                    result += [pp]*(i+1)
    return result

def kronecker(f):
    """Recursive factorization of an univariate polynomial with integer
    coeffs using interpolation.
    """
    def lagrange_base(pos):
        """Compute the base polynomials used for interpolation."""
        l=[]
        for x in pos:
            l.append(Polynomial(coeffs=((S.One, S.One),(-x, S.Zero)),
                                var=f.var, order=f.order))
        b=[]
        for i, x in enumerate(pos):
            p = Polynomial(coeffs=((S.One, S.Zero),), var=f.var, order=f.order)
            for ll in l[:i]+l[i+1:]:
                p *= ll
            c = S.One
            for xx in pos[:i]+pos[i+1:]:
                c *= (x-xx)
            p = Polynomial(coeffs=tuple(map(lambda t:(t[0]/c,)+t[1:], p.coeffs)),
                           var=p.var, order=p.order)
            b.append(p)
        return b

    def combine(divs):
        # Don't try negative divisors for first value.
        lst = map(lambda el: [el], divs[0])
        for choices in divs[1:]:
            lst2 = []
            for el in lst:
                for new in choices:
                    # Also try negative divisors:
                    lst2.append(el + [new])
                    lst2.append(el + [-new])
            lst = lst2
        return lst

    # Half the degree for a possible polynomial divisor g.
    deg = int(f.coeffs[0][1])/2
    # Points for interpolation
    pos = map(Rational, range(0-deg/2, deg+1-deg/2))
    # Reusable Lagrange base polynomials.
    base = lagrange_base(pos)
    # Evaluate at points.
    values = map(f, pos)
    # Look for roots, that is, zeros in values:
    # TODO: Move out of kronecker()!
    lfs = []
    for x, y in zip(pos, values):
        if y == 0:
            lfs.append(Polynomial(coeffs=((S.One, S.One), (-x, S.Zero)),
                                  var=f.var, order=f.order))
    if len(lfs) > 0:
        ff = Polynomial(S.One, ((S.One, S.Zero),), f.var, f.order)
        for lf in lfs:
            ff *= lf
        return lfs + kronecker(div_.mv_int(f, ff)[0][0])
    # All divisors of the values give possible values for g.
    divs = map(integer_divisors, map(int, values))
    # Assemble all possible divisor combination
    combs = combine(divs)
    # Construct candidates for g.
    cands = []
    for comb in combs:
        cand = Polynomial(S.Zero, var=f.var, order=f.order)
        for c,b in zip(comb, base):
            cand += Polynomial(coeffs=tuple([(c*term[0],) + term[1:]
                                             for term in b.coeffs]),
                               var=b.var, order=b.order)
            
        # Filter out constant and non-integer polynomials!
        if not (len(cand.coeffs) == 1
                and cand.coeffs[0][1] is S.Zero):
            if all(map(lambda t:t[0].is_integer, cand.coeffs)):
                cands.append(cand)

    # Make leading coefficient positive:
    for cand in cands:
        if cand.coeffs[0][0] < S.Zero:
            cand = Polynomial(coeffs=tuple([(t[0]*S.NegativeOne,) + t[1:]
                                            for t in cand.coeffs]),
                              var=cand.var, order=cand.order)

    # Filter double entries:
    cands2 = []
    for cand in cands:
        if not cand in cands2:
            cands2.append(cand)
    cands = cands2

    # TODO: Too many candidates?
    # TODO: Use iterators instead of lists!

    for g in cands:
        q, r =  div_.mv_int(f,g)
        if r.sympy_expr is S.Zero:
            return kronecker(q[0]) + kronecker(g)
    else:
        # No divisor found, f irreducible.
        # TODO: Try again with smaller degree divisors?
        return [f]

def kronecker_mv(f):
    """Multivariate polynomial factorization.

    Depends on univariate factorization, and hence is restricted to
    rational coeffs. Expects an instance of Polynomial with at
    least 2 var.

    """
    
    def factor_combinations(lisp, m):
        def recursion(fa, lisp, m):
            if m == 0:
                 yield fa
            else:
                for i, fa2 in enumerate(lisp[0 : len(lisp) + 1 - m]):
                    for el in recursion(fa2 * fa, lisp[i + 1:], m - 1):
                        yield el

        for i, fa in enumerate(lisp[0 : len(lisp) + 1 - m]):
            for el in recursion(fa, lisp[i + 1:], m - 1):
                yield el

    # First sort the var by occuring exponents.
    # Then get degree bound, that is larger than all individual degrees.
    max_exp = {}
    for v in f.var:
        max_exp[v] = 0
    for term in f.coeffs:
        for v, exponent in zip(f.var, term[1:]):
            if exponent > max_exp[v]:
                max_exp[v] = exponent
    new_var = list(f.var)
    new_var.sort(key=lambda v:max_exp[v], reverse=True)
    f = Polynomial(f.sympy_expr, var=new_var, order=f.order)
    d = int(max_exp[f.var[0]]) + 1

    # Now reduce the polynomial f to g in just variable, by the
    # substitution x_i -> y**(d**i)
    g = f.sympy_expr
    y = Symbol('y', dummy=True)
    for i, v in enumerate(f.var):
        g = g.subs(v, y**(d**i))
    g = Polynomial(g, var=y, order=f.order)

    # We can now call the univariate factorization algorithm for g.
    g_factors = uv_int(g)

    # Trial division with all combinations of factors of g.
    tested = []
    result = []
    for m in range(1, len(g_factors)):
        for cand in factor_combinations(g_factors, m):
            if cand in tested:
                continue
            # Inverse reduction
            ff = S.Zero
            for term in cand.coeffs:
                ff_term = term[0]
                y_deg = term[1]
                for v in f.var:
                    v_deg = int(y_deg) % d
                    y_deg = (y_deg - v_deg)/d
                    ff_term *= v**v_deg
                ff += ff_term
            if ff is S.One:
                continue
            candidate = Polynomial(ff, var=f.var, order=f.order)
            q, r = div_.mv_int(f, candidate)
            if r.sympy_expr is S.Zero: # found a factor
                result.append(candidate)
                f = q[0]
            else:
                tested.append(cand)
    if f.sympy_expr is not S.One:
        result.append(f)
    return result

## def mv_int(f):
##     """Complete multivariate polynomial factorization over the integers.

##     Depends on univariate factorization.
##     """

##     var = f.var

##     # TODO: Choose best main variable, such that leading coefficient is 1
##     # or the degree is small!
##     x = var[0]
##     f.var = [x]
##     # Get coeffs, now polynomials in the rest of the var
##     c = map(lambda t:t[0], f.cl)
##     c = map(lambda t: Polynomial(t, var=var[1:]), c)

##     # Get content removed from the polynomial, can be factored seperately
##     content = reduce(gcd_.mv, c)
##     for term, coeff in zip(f.cl, c):
##         q, r = div_.mv(coeff, content)
##         term[0] = q[0].basic
##     f.cl = f.cl # To clear f.basic

##     # Do square-free factorization with main variable
##     a = sqf_.uv_int(f)

##     # Run Wang's algorithm for all factors
##     for i, p in enumerate(a):
##         if p.cl[0][1:] != [0]*len(var): # Filter out constants
##             for pp in wang(p, var):
##                 result += [pp]*(i+1)

##     # Combine result with content's factorization
##     return  mv_int(content) + result

## def wang(f, var):
##     """Following
##     Factoring Multivariate Polynomials Over the Integers
##     Paul S. Wang; Linda Preiss Rothschild
##     Mathematics of Computation, Vol. 29, No. 131. (July, 1975), pp. 935-950

##     """

##     def change(subs):
##         i = random.randint(0, len(subs) + 1)
##         subs[i] += 1
##         return subs

##     # First substitute the var in var with integers,
##     # such that the resulting polynomial is of same degree
##     # and still square free
##     # TODO: Which integers to use?
##     substitutes = [S.Zero]*len(var)
##     while True:
##         lead_coeff = f.cl[0][0]
##         basic = f.basic
##         for i, v in enumerate(var):
##             lead_coeff = lead_coeff.subs(v, substitutes[i])
##         if lead_coeff == S.Zero:
##             substitutes = change(substitutes)
##             continue
##         for i, v in enumerate(var):
##             basic = basic.subs(v, substitutes[i])
##         if basic != sqf_part(basic):
##             substitutes = change(substitutes)
##             continue
##         else: # Substitutes are fine, go on.
##             break

##     # Now factorize the remaining polynomial in one variable.
##     ff = Polynomial(basic, f.var, f.order, f.coeff)
##     factors = uv_int(ff)

##     # When univariate polynomial is irreducible, the original was, too.
##     if len(factors) == 1:
##         return Polynomial(f.basic, f.var + var, f.order, f.coeff)

##     # Now, try to reconstruct multivariate factor candidates from
##     # the univariate ones.
    
