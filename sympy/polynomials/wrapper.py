"""Module providing a user-friendly interface to the polynomial algorithms."""

from sympy.polynomials.base import *
from sympy.polynomials import div_
from sympy.polynomials import factor_
from sympy.polynomials import gcd_
from sympy.polynomials import groebner_
from sympy.polynomials import lcm_
from sympy.polynomials import roots_
from sympy.polynomials import sqf_

def div(f, g, var=None, order=None, coeff=None):
    """Polynomial division of f by g, returns quotients and remainder.

    Univariate and multivariate polynomials are possible. The argument
    g can also be a list of polynomials to be used as divisors. This
    algorithm doesn't stop when the leading terms don't divide, but
    instead tries to reduce smaller terms after that. When dealing
    with multiple var, the monomial ordering in use has
    influence over the result. Optionally, a ring of coefficients can
    be indicated, to restrict computations to this domain.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> div(x**2+6*x+1, 3*x-1)
    (19/9 + (1/3)*x, 28/9)
    >>> div(x**2+6*x+1, 3*x-1, coeff='int')
    (2, 3 + x**2)
    >>> div(2*x**3*y**2 - x*y + y**3, [x-y, y**2], [x,y], 'lex')
    ([-y + 2*y**4 + 2*x*y**3 + 2*x**2*y**2, (-1) + y + 2*y**3], 0)

    """
    f = sympify(f)
    if not isinstance(g, list):
        g = [g]
    g = map(lambda x: sympify(x), g)
    if isinstance(var, Symbol):
        var = [var]
    if var is None:
        var = merge_var(list(f.atoms(type=Symbol)),
                        *[list(g_i.atoms(type=Symbol)) for g_i in g])
    f = Polynomial(f, var=var, order=order)
    g = map(lambda x: Polynomial(x, var=var, order=order), g)
    if coeff is not None:
        if not coeff in coeff_rings:
            raise PolynomialException(
                "%s is not an implemented coefficient ring." % coeff)
        elif coeff == 'int':
            # TODO: Check if given polynomials have integer coeffs?
            q, r = div_.mv_int(f, g)
        else:
            q, r = div_.mv(f, g)
    else: # coeff is None
        q, r = div_.mv(f, g)
    if len(q) == 1:
        return q[0].sympy_expr, r.sympy_expr
    else:
        return map(lambda x: x.sympy_expr, q), r.sympy_expr

def factor(f):
    """Factors the polynomial f over the rationals.

    Example:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> factor(x**4 - 1)
    -(1 + x)*(1 + x**2)*(1 - x)
    >>> factor(x**2 - y**2)
    (y - x)*(-x - y)
    """
    l, f = Polynomial(f).as_integer()
    c, f = f.as_primitive()

    # Get an re-assemble factors:
    result = S.One
    if len(f.var) == 1:
        for p in factor_.uv_int(f):
            result *= p.sympy_expr.expand()
    if len(f.var) > 1:
        for p in factor_.kronecker_mv(f):
            result *= p.sympy_expr.expand()
    return result*c/l

def gcd(f, g, var=None, order=None, coeff=None):
    """Greatest common divisor of two polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> gcd(4*x**2*y, 6*x*y**2, coeff='rat')
    x*y
    >>> gcd(4*x**2*y, 6*x*y**2)
    2*x*y

    """
    f = sympify(f)
    g = sympify(g)
    if var is None:
        var = merge_var(f.atoms(type=Symbol), g.atoms(type=Symbol))
    if isinstance(var, Symbol):
        var = [var]
    if coeff is None:
        atoms = f.atoms().union(g.atoms())
        atoms = filter(lambda a: not a in var, atoms)
        coeff = coeff_ring(atoms)
    if len(var) == 0:
        if coeff == 'int':
            assert isinstance(f, Rational) and isinstance(g, Rational)
            assert f.is_integer and g.is_integer
            return abs(numbers.gcd(f.p, g.p))
        else:
            return S.One
    elif len(var) == 1:
        if coeff == 'int':
            cf, f = Polynomial(f, var=var, order=order).as_primitive()
            cg, g = Polynomial(g, var=var, order=order).as_primitive()
            return abs(numbers.gcd(cf.p, cg.p)) * \
                   gcd_.uv(f, g).sympy_expr
        else:
            r =  gcd_.uv(Polynomial(f, None, var, order),
                         Polynomial(g, None, var, order))
            return r.sympy_expr
    else:
        if coeff == 'int':
            cf, f = Polynomial(f, var=var, order=order).as_primitive()
            cg, g = Polynomial(g, var=var, order=order).as_primitive()
            return abs(numbers.gcd(cf.p, cg.p)) * \
                   gcd_.mv(f, g).sympy_expr
        else:
            r =  gcd_.mv(Polynomial(f, None, var, order),
                         Polynomial(g, None, var, order))
            return r.sympy_expr

def groebner(f, var=None, order=None, reduced=True):
    """Computes the (reduced) Groebner base of given polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> groebner([x**2 + y**3, y**2-x], order='lex')
    [x - y**2, y**3 + y**4]

    """
    if isinstance(f, Basic):
        f = [f]
    f = map(lambda p: sympify(p).expand(), f)
    # filter trivial or double entries
    ff = filter(lambda x: x!=0, f)
    if not ff: # Zero Ideal
        return [S.Zero]
    f = []
    for p in ff:
        if not p in f:
            f.append(p)
    if var is None:
        var = merge_var(*map(lambda p: p.atoms(type=Symbol),f))
    if isinstance(var, Symbol):
        var = [var]
    f = map(lambda p: Polynomial(p, var=var, order=order), f)
    g = groebner_.groebner(f, reduced)
    return map(lambda p: p.sympy_expr, g)

def lcm(f, g, var=None, order=None, coeff=None):
    """Least common divisor of two polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> lcm(4*x**2*y, 6*x*y**2, coeff='rat')
    x**2*y**2
    >>> lcm(4*x**2*y, 6*x*y**2)
    12*x**2*y**2
    """
    f = sympify(f)
    g = sympify(g)
    if var is None:
        var = merge_var(f.atoms(type=Symbol), g.atoms(type=Symbol))
    if isinstance(var, Symbol):
        var = [var]
    if coeff is None:
        atoms = f.atoms().union(g.atoms())
        atoms = filter(lambda a: not a in var, atoms)
        coeff = coeff_ring(atoms)
    if len(var) == 0:
        if coeff == 'int':
            assert isinstance(f, Rational) and isinstance(g, Rational)
            assert f.is_integer and g.is_integer
            if f == S.Zero or g == S.Zero:
                return S.Zero
            return abs(f*g / numbers.gcd(f.p, g.p))
        else:
            return S.One
    if len(var) == 1:
        if coeff == 'int':
            cf, f = Polynomial(f, var=var, order=order).as_primitive()
            cg, g = Polynomial(g, var=var, order=order).as_primitive()
            if cf is S.Zero or cg is S.Zero:
                return S.Zero
            return abs(cf*cg / numbers.gcd(cf.p, cg.p)) * \
                   lcm_.uv(f, g).sympy_expr
        else:
            r = lcm_.uv(Polynomial(f, var=var, order=order),
                        Polynomial(g, var=var, order=order))
            return r.sympy_expr
    else:
        if coeff == 'int':
            cf, f = Polynomial(f, var=var, order=order).as_primitive()
            cg, g = Polynomial(g, var=var, order=order).as_primitive()
            if cf is S.Zero or cg is S.Zero:
                return S.Zero
            return abs(cf*cg / numbers.gcd(cf.p, cg.p)) * \
                   lcm_.mv(f, g).sympy_expr
        else:
            r = lcm_.mv(Polynomial(f, var=var, order=order),
                        Polynomial(g, var=var, order=order))
            return r.sympy_expr

def real_roots(f, a=None, b=None):
    """Returns the number of unique real roots of f in the interval (a, b].

    Assumes an univariate, square-free polynomial with real coeffs.
    The boundaries a and b can be omitted to check the whole real line.

    Examples:
    >>> x = Symbol('x')
    >>> real_roots(x**2 - 1)
    2
    >>> real_roots(x**2 - 1, 0, 2)
    1

    """
    f = sympify(f)
    if a is not None:
        a = sympify(a)
    if b is not None:
        b = sympify(b)
    var = list(f.atoms(type=Symbol))
    if len(var) == 1:
        var = var[0]
    else:
        raise PolynomialException('Not an univariate polynomial.')
    ss = roots_.sturm(Polynomial(f))
    return roots_.real_roots(ss, a, b)

def resultant(f, g, x, method='bezout'):
    """Computes resultant of two polynomials in one variable. This
       method can be used to verify if selected polynomials have
       common root withot factoring them or computing any GCD's.

       It can be also be used for variable elemination in case of
       bivariate polynomials. Just assume that one of the var
       is a parameter and compute resultant with respect to the other
       one and you will get univariate polynomial in single variable.

       For now two algorithm have been implemented, based on
       Sylvester and Bezout matrices. The default is Bezout.

       TODO: Make Bezout approach run in O(s**2). Currently
             it is O(s**3), where s = max(deg(f), deg(g)).
    """

    from sympy.matrices import zero

    fp = Polynomial(f, var=x).coeffs
    gp = Polynomial(g, var=x).coeffs

    m, n = int(fp[0][1]), int(gp[0][1])

    if method is 'sylvester':
        M = zero(m+n)

        for i in range(n):
            for coeff, j in fp:
                M[i, m-int(j)+i] = coeff

        for i in range(m):
            for coeff, j in gp:
                M[i+n, n-int(j)+i] = coeff

        return M.det()
    elif method is 'bezout':
        if n > m:
            s, fp, gp = n, gp, fp
        else:
            s = m

        p = [Basic.Zero()] * (s+1)
        q = [Basic.Zero()] * (s+1)

        for coeff, j in fp:
            p[int(j)] = coeff

        for coeff, j in gp:
            q[int(j)] = coeff

        M = zero(s)

        for i in range(s):
            for j in range(i, s):
                z = 1 + min(i, s-1-j)
                terms = [Basic.Zero()] * z

                for k in range(z):
                    terms[k] = p[j+k+1]*q[i-k] - p[i-k]*q[j+k+1]

                M[i, j] = M[j, i] = Add(*terms)

        if (s*(s-1)/2) % 2 == 0:
            b = p[-1]**abs(n-m)
        else:
            b = -p[-1]**abs(n-m)

        return (1/b) * M.det()
    else:
        raise ValueError("Invalid method: '%s'" % method)

def roots(f, var=None):
    """Compute the roots of an univariate polynomial.

    The coeff argument determines the type of the roots to look
    for. The supported types include 'int', 'rat', 'real' and 'cplx'
    and the coeffs of given polynomials are assumed to be of
    this type.

    Examples:
    >>> x = Symbol('x')
    >>> roots(x**2 - 1)
    [-1, 1]
    
    """

    f = Polynomial(f, var=var)
    if len(f.var) == 0:
        return []
    if len(f.var) > 1:
        raise PolynomialException(
            'Multivariate polynomials not supported.')
    # Determine type of coeffs (for factorization purposes)
    symbols = f.sympy_expr.atoms(type=Symbol)
    symbols = filter(lambda a: not a in f.var, symbols)
    if symbols:
        coeff = 'sym'
    else:
        coeff = coeff_ring(get_numbers(f.sympy_expr))
        
    if coeff == 'rat':
        l, f = f.as_integer()
        coeff = 'int'
    if coeff in ['int', 'rat']:
        c, f = f.as_primitive()
        # Hacks for special cases (where factorization would do harm)
        # TODO: Support these cases differently, in roots_.uv
        if len(f.coeffs) == 1:
            if f.coeffs[0][1] == 0:
                return []
            else:
                return[S.Zero]
        if len(f.coeffs) == 2:
            result = []
            if f.coeffs[1][1] > 0:
                result.append(S.Zero)
                f = Polynomial(coeffs=tuple([(t[0], t[1] - f.coeffs[1][1])
                                             for t in f.coeffs]),
                               var=f.var, order=f.order)
            return result + roots_.uv(f)    
        else:
            # Factor without multiplicity.
            factors = []
            factors = filter(lambda p: not p in factors,
                             factor_.uv_int(f))
    else:
        factors = [f]
        
    # Now check for roots in each factor
    result = []
    for p in factors:
        result += roots_.uv(p)
        
    return result

def solve_system(eqs, var=None):
    """Solve a system of polynomial equations.

    Only works for finite sets of solutions, which is not
    tested. Currently, some equations can't be solved.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> f = y - x           
    >>> g = x**2 + y**2 - 1 
    >>> solve_system([f, g])
    [[-1/2*2**(1/2), -1/2*2**(1/2)], [(1/2)*2**(1/2), (1/2)*2**(1/2)]]
    
    """
    if not isinstance(eqs, list):
        eqs = [eqs]
    if var is None:
        var = merge_var(*[eq.atoms(type=Symbol) for eq in eqs])
    eqs = map(lambda eq:Polynomial(eq, None, var, 'lex'), eqs)
    return roots_.equation_system(eqs)

def sqf(f, var=None):
    """Computes the square-free decomposition of 'f'.

    Only works for univariate polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> sqf(2*x**3 + 2*x**2)
    x**2*(2 + 2*x)

    """
    a = sqf_.uv(Polynomial(f, var=var))

    result = 1
    for i, p in enumerate(a):
        result *= (p.sympy_expr)**(i+1)
    return result

def sqf_part(f, var=None):
    """Computes the square-free part of f.

    Only works for univariate polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> sqf_part(2*x**3 + 2*x**2)
    2*x + 2*x**2

    """
    return (sqf_.uv_part(Polynomial(f, var=var))).sympy_expr
