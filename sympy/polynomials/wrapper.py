"""Module providing a user-friendly interface to the polynomial algorithms."""

from sympy.polynomials.base import *
from sympy.polynomials import div_
from sympy.polynomials import factor_
from sympy.polynomials import groebner_
from sympy.polynomials import roots_

def div(f, g, var=None, order=None, coeff=None):
    """Polynomial division of f by g, returns quotients and remainder.

    Univariate and multivariate polynomials are possible. The argument
    g can also be a list of polynomials to be used as divisors. This
    algorithm doesn't stop when the leading terms don't divide, but
    instead tries to reduce smaller terms after that. When dealing
    with multiple var, the monomial ordering in use has
    influence over the result. Optionally, a ring of coefficients can
    be indicated, to restrict computations to this domain (the integers.)

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
    q, r = div_.div(f, g, var, order, coeff)
    
    if isinstance(q, list):
        return map(lambda x: x.sympy_expr, q), r.sympy_expr
    else:
        return q.sympy_expr, r.sympy_expr

def factor(f, var=None, order=None):
    """Factors the polynomial f over the rationals.

    Example:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> factor(x**4 - 1)
    -(1 + x)*(1 + x**2)*(1 - x)
    >>> factor(x**2 - y**2)
    (y - x)*(-x - y)
    """
    # Re-assemble factors:
    return Mul(*[ff.sympy_expr for ff in factor_.factor(f, var, order)])

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

    # First check if the polynomials have integer coefficients.
    if coeff is None:
        coeff_f = coeff_ring(get_numbers(f))
        coeff_g = coeff_ring(get_numbers(g))
        if coeff_f == 'int' and coeff_g == 'int':
            coeff = 'int'
        if var is not None:
            if isinstance(var, Symbol):
                var = [var]
            if [v for v in f.atoms(type=Symbol).union(g.atoms(type=Symbol))
                if v not in var]:
                coeff = 'sym'
    
    return div_.gcd(f, g, var, order, coeff).sympy_expr

def groebner(f, var=None, order=None, reduced=True):
    """Computes the (reduced) Groebner base of given polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> groebner([x**2 + y**3, y**2-x], order='lex')
    [x - y**2, y**3 + y**4]

    """
    g = groebner_.groebner(f, var, order, reduced)
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
    # First check if the polynomials have integer coefficients.
    if coeff is None:
        coeff_f = coeff_ring(get_numbers(f))
        coeff_g = coeff_ring(get_numbers(g))
        if coeff_f == 'int' and coeff_g == 'int':
            coeff = 'int'
        if var is not None:
            if isinstance(var, Symbol):
                var = [var]
            if [v for v in f.atoms(type=Symbol).union(g.atoms(type=Symbol))
                if v not in var]:
                coeff = 'sym'
    
    return div_.lcm(f, g, var, order, coeff).sympy_expr

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

def sqf(f, var=None, order=None, coeff=None):
    """Computes the square-free decomposition of 'f'.

    Only works for univariate polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> sqf(2*x**3 + 2*x**2)
    x**2*(2 + 2*x)

    """
    a = factor_.sqf(f, var, order, coeff)
    result = S.One
    for i, p in enumerate(a):
        result *= (p.sympy_expr)**(i+1)
    return result

def sqf_part(f, var=None, order=None):
    """Computes the square-free part of f.

    Only works for univariate polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> sqf_part(2*x**3 + 2*x**2)
    2*x + 2*x**2

    """
    return factor_.sqf_part(f, var=var, order=order).sympy_expr
