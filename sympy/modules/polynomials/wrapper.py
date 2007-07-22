"""Module providing a user-friendly interface to the polynomial algorithms."""

#from sympy.core.functions import diff
#import sympy.modules.matrices

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import div_
from sympy.modules.polynomials import factor_
from sympy.modules.polynomials import gcd_
from sympy.modules.polynomials import groebner_
from sympy.modules.polynomials import lcm_
from sympy.modules.polynomials import roots_
from sympy.modules.polynomials import sqf_

# TODO: Remove/incorporate in Polynomial
def coeff(poly, x, n):
    """Returns the coefficient of x**n in the polynomial"""
    assert ispoly(poly,x)
    p = Polynomial(poly, x)
    t = filter(lambda t:t[1]==n,p.cl)
    if len(t) == 0:
        return S.Zero
    else:
        return t[0][0]

def div(f, g, var=None, order=None, coeff=None):
    """Polynomial division of f by g, returns quotients and remainder.

    Univariate and multivariate polynomials are possible. The argument
    g can also be a list of polynomials to be used as divisors. This
    algorithm doesn't stop when the leading terms don't divide, but
    instead tries to reduce smaller terms after that. When dealing
    with multiple variables, the monomial ordering in use has
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
    f = Basic.sympify(f)
    if not isinstance(g, list):
        g = [g]
    g = map(lambda x: Basic.sympify(x), g)
    if not isinstance(var, list):
        var = [var]
    if len(var) > 0 and var[0] is None:
        var = merge_var(list(f.atoms(type=Symbol)),
                        *[list(g_i.atoms(type=Symbol)) for g_i in g])
    f = Polynomial(f, var, order, coeff)
    g = map(lambda x: Polynomial(x, var, order, coeff), g)
    if coeff is not None:
        if not coeff in coeff_rings:
            raise PolynomialException(
                "%s is not an implemented coefficient ring." % coeff)
        elif coeff == 'int':
            # TODO: Check if given polynomials have integer coeffs?
            q, r = div_.mv_int(f, g)
    else: # coeff is None
        q, r = div_.mv(f, g)
    if len(q) == 1:
        return q[0].basic, r.basic
    else:
        return map(lambda x: x.basic, q), r.basic

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
    f = Polynomial(f, coeff='int')

    # Compute lcm of denominators in coefficients:
    l = 1
    for term in f.cl:
        if not isinstance(term[0], Rational):
            raise PolynomialException('Non-rational coefficient!')
        l = l*term[0].q / numbers.gcd(l, term[0].q)

    # Make f a polynomial in integer coefficients:
    f.cl = map(lambda t:[t[0]*Rational(l)] + t[1:], f.cl)

    # Remove content before factorization:
    c = f.content()
    f.cl = map(lambda t:[t[0]/c] + t[1:], f.cl)

    # Get an re-assemble factors:
    result = S.One
    if len(f.var) == 1:
        for p in factor_.uv_int(f):
            result *= p.basic.expand()
    if len(f.var) > 1:
        for p in factor_.kronecker_mv(f):
            result *= p.basic.expand()
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
    f = Basic.sympify(f)
    g = Basic.sympify(g)
    if isinstance(var, Symbol):
        var = [var]
    if var is None:
        var = merge_var(list(f.atoms(type=Symbol)), list(g.atoms(type=Symbol)))
    if coeff is None:
        atoms = list(f.atoms())
        atoms += filter(lambda a: not a in atoms, list(g.atoms()))
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
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            return abs(numbers.gcd(cf.p, cg.p)) * \
                   gcd_.uv(f, g).basic
        else:
            r =  gcd_.uv(Polynomial(f, var, order, coeff),
                         Polynomial(g, var, order, coeff))
            return r.basic
    else:
        if coeff == 'int':
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            return abs(numbers.gcd(cf.p, cg.p)) * \
                   gcd_.mv(f, g).basic
        else:
            r =  gcd_.mv(Polynomial(f, var, order, coeff),
                         Polynomial(g, var, order, coeff))
            return r.basic

def groebner(f, var=None, order=None, coeff=None, reduced=True):
    """Computes the (reduced) Groebner base of given polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> groebner([x**2 + y**3, y**2-x], order='lex')
    [x - y**2, y**3 + y**4]

    """
    if isinstance(f, Basic):
        f = [f]
    f = map(lambda p: Basic.sympify(p).expand(), f)
    # filter trivial or double entries
    ff = filter(lambda x: x!=0, f)
    if not ff: # Zero Ideal
        return [S.Zero]
    f = []
    for p in ff:
        if not p in f:
            f.append(p)
    if isinstance(var, Symbol):
        var = [var]
    if var is None:
        var = merge_var(*map(lambda p: list(p.atoms(type=Symbol)),f))
    f = map(lambda p: Polynomial(p, var, order, coeff), f)
    g = groebner_.groebner(f, reduced)
    return map(lambda p: p.basic, g)

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
    f = Basic.sympify(f)
    g = Basic.sympify(g)
    if isinstance(var, Symbol):
        var = [var]
    if var is None:
        var = merge_var(list(f.atoms(type=Symbol)), list(g.atoms(type=Symbol)))
    if coeff is None:
        atoms = list(f.atoms())
        atoms += filter(lambda a: not a in atoms, list(g.atoms()))
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
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            if cf == S.Zero or cg == S.Zero:
                return S.Zero
            return abs(cf*cg / numbers.gcd(cf.p, cg.p)) * \
                   lcm_.uv(f, g).basic
        else:
            r = lcm_.uv(Polynomial(f, var, order, coeff),
                        Polynomial(g, var, order, coeff))
            return r.basic
    else:
        if coeff == 'int':
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            if cf == S.Zero or cg == S.Zero:
                return S.Zero
            return abs(cf*cg / numbers.gcd(cf.p, cg.p)) * \
                   lcm_.mv(f, g).basic
        else:
            r = lcm_.mv(Polynomial(f, var, order, coeff),
                        Polynomial(g, var, order, coeff))
            return r.basic

def real_roots(f, a=None, b=None):
    """Returns the number of unique real roots of f in the interval (a, b].

    Assumes an univariate, square-free polynomial with real coefficients.
    The boundaries a and b can be omitted to check the whole real line.

    Examples:
    >>> x = Symbol('x')
    >>> real_roots(x**2 - 1)
    2
    >>> real_roots(x**2 - 1, 0, 2)
    1

    """
    f = Basic.sympify(f)
    if a is not None:
        a = Basic.sympify(a)
    if b is not None:
        b = Basic.sympify(b)
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
       bivariate polynomials. Just assume that one of the variables
       is a parameter and compute resultant with respect to the other
       one and you will get univariate polynomial in single variable.

       For now two algorithm have been implemented, based on
       Sylvester and Bezout matrices. The default is Bezout.

       TODO: Make Bezout approach run in O(s**2). Currently
             it is O(s**3), where s = max(deg(f), deg(g)).
    """

    from sympy.modules.matrices import zero

    fp = coeff_list(f, x)
    gp = coeff_list(g, x)

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

def roots(f, var=None, coeff=None, verbose=False):
    """Compute the roots of an univariate polynomial.

    The coeff argument determines the type of the roots to look
    for. The supported types include 'int', 'rat', 'real' and 'cplx'
    and the coefficients of given polynomials are assumed to be of
    this type.

    Examples:
    >>> x = Symbol('x')
    >>> roots(x**2 - 1)
    [-1, 1]
    
    """

    f = Polynomial(f, var)
    if len(f.var) != 1:
        raise PolynomialException(
            'Multivariate polynomials not supported.')
    if coeff is None:
        # Allow all roots by default.
        coeff = coeff_rings[-1]
    # Determine type of coefficients (for factorization purposes)
    symbols = list(f.basic.atoms(type=Symbol))
    symbols = filter(lambda a: not a in f.var, symbols)
    if symbols:
        coeff2 = 'sym'
    else:
        coeff2 = coeff_ring(get_numbers(f.basic))
    f.coeff = coeff2
    
    if f.coeff == 'rat':
        # Compute lcm of denominators in coefficients:
        l = 1
        for term in f.cl:
            if not isinstance(term[0], Rational):
                raise PolynomialException('Non-rational coefficient!')
            l = l*term[0].q / numbers.gcd(l, term[0].q)

        # Make f a polynomial in integer coefficients:
        f.cl = map(lambda t:[t[0]*Rational(l)] + t[1:], f.cl)
        f.coeff = 'int'
    if f.coeff in ['int', 'rat']:
        # Remove content before factorization:
        c = f.content()
        f.cl = map(lambda t:[t[0]/c] + t[1:], f.cl)
        # Hacks for special cases (where factorization would do harm)
        # TODO: Support these cases differently, in roots_.uv
        if len(f.cl) == 1:
            if f.cl[0][1] == 0:
                return []
            else:
                return[S.Zero]
        if len(f.cl) == 2:
            f.coeff = coeff
            result = []
            if f.cl[1][1] > 0:
                result.append(S.Zero)
                f.cl = map(lambda t:[t[0], t[1]-f.cl[1][1]], f.cl)
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
        p.coeff = coeff
        result += roots_.uv(p)
        
    return result

def solve_system(eqs, var=None, coeff=None):
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
        var = merge_var(*[list(eq.atoms(type=Symbol)) for eq in eqs])
    if coeff is None:
        # Allow all roots by default.
        coeff = coeff_rings[-1]
    eqs = map(lambda eq:Polynomial(eq, var=var, order='lex', coeff=coeff), eqs)
    return roots_.equation_system(eqs)

def sqf(f, var=None):
    """Computes the square-free decomposition of 'f'.

    Only works for univariate polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> sqf(2*x**3 + 2*x**2)
    x**2*(2 + 2*x)

    """
    f = Basic.sympify(f)
    if isinstance(var, Symbol):
        var = [var]
    if var is None:
        var = list(f.atoms(type=Symbol))
    if len(var) != 1:
        raise PolynomialException('Not an univariate polynomial.')

    a = sqf_.uv(Polynomial(f, var))

    result = 1
    for i, p in enumerate(a):
        result *= (p.basic)**(i+1)
    return result

def sqf_part(f, var=None):
    """Computes the square-free part of f.

    Only works for univariate polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> sqf_part(2*x**3 + 2*x**2)
    2*x + 2*x**2

    """
    f = Basic.sympify(f)
    if isinstance(var, Symbol):
        var = [var]
    if var is None:
        var = list(f.atoms(type=Symbol))
    if len(var) != 1:
        raise PolynomialException('Not an univariate polynomial.')

    return (sqf_.uv_part(Polynomial(f, var))).basic
