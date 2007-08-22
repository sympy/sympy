"""Various algorithms for the factorization of polynomials."""

from sympy import ntheory

from sympy.polynomials.base import *
from sympy.polynomials import div_, fast

def sqf(f, var=None, order=None, coeff=None):
    """Square-free decomposition.

    Usage:
    ======
        Computes a decomposition of f in a1 * a2**2 * ... * an**n,
        where the ai are pairwise prime and square-free polynomials.

        The input is assumed to be a univariate polynomial, either as
        a SymPy expression or an instance of Polynomial. In the first
        case, you can optionally specify the variables and monomial
        order with the arguments 'var' and 'order'.

        If the argument 'coeff' is set to 'int', the constant factors
        are redistributed to the different ai, so that they have
        integer coefficients. Otherwise, they are made monic.

        A list is returned, with an instance of Polynomial at each
        index, which represents the multiplicity (except beginning
        with 1 instead of 0).
    
    Examples:
    =========
        >>> x = Symbol('x')
        >>> a = sqf(3 - 12*x - 4*x**3 + 4*x**4 + 13*x**2)
        >>> for i, f in enumerate(a): print (i + 1), f
        1 12 + 4*x**2
        2 (-1/2) + x
        >>> b = sqf(3 - 12*x - 4*x**3 + 4*x**4 + 13*x**2, coeff='int')
        >>> for i, f in enumerate(b): print (i + 1), f
        1 3 + x**2
        2 (-1) + 2*x
        
    References:
    ===========
        Gathen, Gerhard: Modern Computer Algebra,
        Cambridge University Press, 1. edition, p. 371

    Also see L{sqf_part}, L{factor}.

    """

    if not isinstance(f, Polynomial):
        f = Polynomial(f, var=var, order=order)

    # Check for constant polynomials:
    if f.var == [] or f.coeffs[0][1] is S.Zero:
        return [f]
    
    f = [f]
    while f[-1].coeffs[0][1] is not S.Zero:
        f.append(div_.gcd(f[-1], f[-1].diff(f[-1].var[0])))
    g = []
    for i in range(1, len(f)):
        g.append(div_.div(f[i-1], f[i])[0])
    a = []
    for i in range(0, len(g)-1):
        a.append(div_.div(g[i], g[i+1])[0])
    a.append(g[-1])
    
    if coeff == 'int': # Redistribute the constants.
        ca = int(a[0].content())
        c = ca
        for i, p in enumerate(a[1:]):
            # Compute lcm of denominators in coeffs:
            l, a[i+1] = p.as_integer()
            c /= int(l)**(i+2)
        ca = Integer(ca/c)
        a[0] = Polynomial(coeffs=tuple([(t[0]/ca,) + t[1:]
                                        for t in a[0].coeffs]),
                          var=a[0].var, order=a[0].order)
    return a


def sqf_part(f, var=None, order=None):
    """Returns the square-free part of f.

    Usage:
    ======
        Computes the square-free part of f.

        The input is assumed to be a univariate polynomial, either as
        a SymPy expression or an instance of Polynomial. In the first
        case, you can optionally specify the variables and monomial
        order with the arguments 'var' and 'order'.

        The result is returned as an instance of Polynomial.
    
    Examples:
    =========
        >>> x = Symbol('x')
        >>> print sqf_part(2*x**3 + 2*x**2)
        2*x + 2*x**2

    References:
    ===========
        Gathen, Gerhard: Modern Computer Algebra,
        Cambridge University Press, 1. edition, p. 370

    Also see L{sqf}.

    """

    if not isinstance(f, Polynomial):
        f = Polynomial(f, var=var, order=order)

    # The gcd of f with its derivative is the multiple part.
    ff = div_.gcd(f, f.diff(f.var[0]))
    return div_.div(f, ff)[0]


def factor(f, var=None, order=None):
    """Factorization of polynomials over the rationals.

    Usage:
    ======
        As input, a polynomial is taken, either as a SymPy expression
        or as an instance of Polynomial. In the first case, the
        variables and monomial order can optionally be specified
        through the arguments 'var' and 'order'.

        The result is a list containing all the irreducible factors as
        instances of Polynomial, the first element being a constant
        factor. 

    Notes:
    ======
        In the univariate case, the first step is a square-free
        decomposition. Then, for each factor, we get all the linear
        factors by looking for rational roots. The remainder is
        finally factored by a method due to Kronecker, which uses
        polynomial interpolation to guess divisors. Due to this
        method, the univariate factoring is rather limited, but
        planned to be replaced soon.

        The multivariate polynomials are reduced to univariate ones,
        also due to a method by Kronecker. The univariate factors are
        then re-assembled and reformed to potential divisors of f, for
        trial division.

        Factorization over the integers and rationals is equivalent,
        since you can always multiply with the common denominator and
        then remove the content of the polynomial.

    Examples:
    =========
        # Order of terms is no longer deterministic,
        # implement cmp for Polynomial?
        #>>> x, y = symbols('xy')
        #>>> for f in factor(2*x**4 - 2): print f
        #2
        #(-1) + x
        #1 + x
        #1 + x**2
        #>>> for f in factor(4*x**2/3 - y**2/3): print f
        #1/3
        #y + 2*x
        #-y + 2*x

    References:
    ===========
        Buchberger, Collins, Loos, Albrecht: Computer Algebra,
        ACM SIGSAM Bulletin 1982

    Also see L{sqf}, L{kronecker}, L{kronecker_mv}, L{roots_.rat_roots}.
    
    """

    from sympy.polynomials import roots_
    
    if not isinstance(f, Polynomial):
        f = Polynomial(f, var=var, order=order)

    # Make it a primitive polynomial in integer coefficients.
    denom, f = f.as_integer()
    content, f = f.as_primitive()
    result = [Polynomial(content/denom, var=f.var, order=f.order)]

    if len(f.var) == 0 or f.coeffs[0][1:] == tuple([S.Zero]*len(f.var)):
        return result
    elif len(f.var) == 1:
        factors = fast.intpoly.factor(Polynomial2IntPoly(f))
        result[0] *= IntPoly2Polynomial(factors[0][0], f.var, f.order)
        for ff in factors[1:]:
            result += [IntPoly2Polynomial(ff[0], f.var, f.order)]*ff[1]
    else: # len(f.var) > 1
        factors = kronecker_mv(f)
        result[0] *= factors[0]
        result += factors[1:]
    return result


def kronecker(f):
    """One step in univariate factorization, see L{factor}."""
    
    def lagrange_base(pos):
        """Compute the base polynomials used for Lagrange interpolation.

        They are constructed such that they evaluate to 1 at their
        position but 0 at all other points.

        """

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
        """Combine the divisors of all coefficients."""
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
    # All divisors of the values give possible values for g.
    divs = map(integer_divisors, map(int, values))
    # Assemble all possible divisor combinations.
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
        q, r =  div_.div(f, g, coeff='int')
        if r.sympy_expr is S.Zero:
            return kronecker(q) + kronecker(g)
    else:
        # No divisor found, f irreducible.
        # TODO: Try again with divisors of smaller degree?
        return [f]


def kronecker_mv(f):
    """One step in multivariate factorization, see L{factor}."""

    # Given a list of factors, re-assemble all m-subsets of them.
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

    # First sort the variables by occuring exponents.
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
    g_factors = factor(g)[1:] # Don't use constant factor.
    constant_factor = Polynomial(S.One, var=f.var, order=f.order)

    # Trial division with all combinations of factors of g.
    tested = []
    result = []
    for m in range(1, len(g_factors)/2 + 1):
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
            # Make leading_coefficient positive:
            if candidate.coeffs[0][0] < 0:
                candidate = Polynomial(
                    coeffs=[(-t[0],)+t[1:] for t in candidate.coeffs],
                    var=candidate.var, order=candidate.order)
            q, r = div_.div(f, candidate, coeff='int')
            if r.sympy_expr is S.Zero: # found a factor
                result.append(candidate)
                f = q
            else:
                tested.append(cand)
            # Check if f is constant.
            if f.coeffs[0][1:] == tuple([S.Zero]*len(f.var)):
                constant_factor = f
                f = Polynomial(S.One, var=f.var, order=f.order)
                break
        if f.sympy_expr is S.One:
            break
    if f.sympy_expr is not S.One:
        result.append(f)

    return [constant_factor] + result
