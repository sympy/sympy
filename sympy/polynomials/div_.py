"""Polynomial division algorithms for use with class Polynomial"""

from sympy.polynomials.base import *

def div(f, g, var=None, order=None, coeff=None):
    """Division with remainder.

    Usage:
    ======
        The input consists of a polynomial f and either one or a list
        of polynomials g. When these are just SymPy expressions, you
        can additionally specify the variables and monomial orders
        with 'var' and 'order', respectively. Only f is checked for
        the input types, the rest is assumed to match.

        If 'coeff' is set to 'int', only term divisions with proper
        coefficient divisions are allowed. That is, 3*x divides 6*x*y,
        but not 2*x**2.

        The output's type is Polynomial, but there is a wrapper, see
        L{wrapper.div} With only one element in g given, the resulting
        list of quotients is unpacked, also. The second output is the
        remainder.

    Notes:
    ======
        The loop then iterates over all terms of f, checking if any of
        the elements in g (or just g if it is the sole divisor) have a
        leading term dividing the current term of f. If yes, the
        proper multiple of the element of g is subtracted from f, so
        that the term is eliminated, otherwise added to the remainder,
        until f is 0.

        This way, the algorithm doesn't stop, when the leading terms
        of g don't divide the leading term of f, but rather try to
        reduce all of f's other terms. Of course, the known univariate
        single-polynomial division with remainder is a special case of
        this function.

        The result can depend on the order of the elements in g. For
        unicity, you need that their form a Groebner base of the ideal
        they generate, see L{groebner_.groebner}.

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> q, r = div(x**2 + 6*x + 1, 3*x - 1)
        >>> print q
        19/9 + (1/3)*x
        >>> print r
        28/9
        >>> q, r = div(x**2 + 6*x + 1, 3*x - 1, coeff='int')
        >>> print q
        2
        >>> print r
        3 + x**2
        >>> q, r = div(2*x**3*y**2 - x*y + y**3, [x-y, y**2], [x,y], 'lex')
        >>> print q[0]
        -y + 2*y**4 + 2*x*y**3 + 2*x**2*y**2
        >>> print q[1]
        (-1) + y + 2*y**3
        >>> print r
        0

    References:
    ===========
        Cox, Little, O'Shea: Ideals, Varieties and Algorithms,
        Springer, 2. edition, p. 62

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

    while f.sympy_expr is not S.Zero:
        for g_i in g:
            if g_i.sympy_expr is S.Zero: # Avoid division by 0.
                continue
            # Check if leading term of f is divisible by that of g_i.
            # When coeff equals 'int', check the coefficient's
            # divisibility, too.
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

    Usage:
    ======
        Input are two polynomials, either as SymPy expressions or as
        instances of Polynomial. In the first case, the variables and
        monomial order can be specified with 'var' and 'order',
        respectively.

        If 'coeff' is set to 'int', the content of each polynomial is
        checked and their gcd multiplied to the result. Otherwise it
        is monic, that is, of leading coefficient 1.

        The output's type is Polynomial, but there is a wrapper, see
        L{wrapper.gcd}.

    Notes:
    ======
        With only one variable, euclid's algorithm is used directly,
        which is reasonably fast. But in the multivariate case, we
        have to compute the gcd using the least common multiple, which
        relies on Groebner bases. This is based on the formula:
        f*g = gcd(f, g)*lcm(f, g)

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> print gcd(4*x**2*y, 6*x*y**2)
        x*y
        >>> print gcd(4*x**2*y, 6*x*y**2, coeff='int')
        2*x*y
        
    References:
    ===========
        Cox, Little, O'Shea: Ideals, Varieties and Algorithms,
        Springer, 2. edition, p. 41 & p. 187

    See also L{div}, L{lcm}.

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
        cf, f = f.as_primitive()
        cg, g = g.as_primitive()
        c = Integer(numbers.gcd(int(cf), int(cg)))
    else:
        c = S.One

    if len(f.var) == 0: # Constant result.
        return Polynomial(c, var=f.var, order=f.order)
    elif len(f.var) == 1: # Use euclidean algorithm.
        while g.sympy_expr is not S.Zero:
            # Remove leading coefficient, to simplify computations.
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
    """Least common multiple.

    Usage:
    ======
        Input are two polynomials, either as SymPy expressions or as
        instances of Polynomial. In the first case, the variables and
        monomial order can be specified with 'var' and 'order',
        respectively.

        If 'coeff' is set to 'int', the content of each polynomial is
        checked and their lcm multiplied to the result. Otherwise it
        is monic, that is, of leading coefficient 1.

        The output's type is Polynomial, but there is a wrapper, see
        L{wrapper.lcm}.

    Notes:
    ======
        With only one variable, the gcd is used to get the lcm from
        the product, via f*g = gcd(f, g)*lcm(f, g).

        In the univariate case, we compute the unique generator of the
        intersection of the two ideals, generated by f and g. This is
        done by computing a lexicographic Groebner base of
        [t*f. (t-1)*g], with t a dummy variable at the first place,
        then filtering trough the base elements not containing t.

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> print lcm(4*x**2*y, 6*x*y**2)
        x**2*y**2
        >>> print lcm(4*x**2*y, 6*x*y**2, coeff='int')
        12*x**2*y**2

    References:
    ===========
        Cox, Little, O'Shea: Ideals, Varieties and Algorithms,
        Springer, 2. edition, p. 187

    See also L{div}, L{gcd}.
        
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
