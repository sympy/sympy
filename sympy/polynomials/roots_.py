"""Algorithms to determine the roots of polynomials"""

from sympy import exp, I, pi, sqrt 

from sympy.polynomials.base import *
from sympy.polynomials import div_, groebner_

def cubic(f):
    """Returns the real or complex roots of a cubic polynomial.

    Works for univariate instances of Polynomial only."""
    # Get monic polynomial, for p-q formula
    f = f.as_monic()[1]
    
    a = f.nth_coeff(2)
    b = f.nth_coeff(1)
    c = f.nth_coeff(0)

    # Substitute variable to get depressed cubic: t**3 + p * t + q
    p = b + -a**2/3
    q = c + Rational(2, 27)*a**3 - a*b/9

    if p == 0: # Solve special cases:
        if q == 0:
            return [-a/3]
        else:
            u1 = q**Rational(1, 3)
    else: 
        u1 = (q/2 + sqrt(q**2/4 + p**3/27))**Rational(1, 3)

    u2 = u1*(Rational(-1, 2) + I*sqrt(3)/2)
    u3 = u1*(Rational(-1, 2) - I*sqrt(3)/2)

    return map(lambda u: (p/(u*3) - u - a/3).expand(), [u1, u2, u3])

def equation_system(eqs):
    """Solves a system of polynomial equation by variable elimination.

    Assumes to get a list of instances of Polynomial, with matching
    var. Only works for zero-dimensional varieties, that is, a
    finite number of solutions (untested!). 
    """
    def is_uv(f):
        """Is an instance of Polynomial univariate in its last variable?
        """
        for term in f.coeffs:
            for exponent in term[1:-1]:
                if exponent > 0:
                    return False
        return True

    # First compute a Groebner base with the polynomials,
    # with lexicographic ordering, so that the last polynomial is
    # univariate and can be solved.
    for i, f in enumerate(eqs):
        eqs[i] = Polynomial(f.sympy_expr, var=f.var, order='lex')
    gb = groebner_.groebner(eqs)

    # Now filter the the base elements, to get only the univariate
    eliminated = filter(is_uv, gb)
    if len(eliminated) != 1:
        raise PolynomialException("System currently not solvable.")

    # Try to solve the polynomials with var eliminated.
    # TODO: Use another function in roots_ (to be written)
    f = eliminated[0]
    from sympy.polynomials import roots
    partial_solutions = roots(f.sympy_expr, f.var[-1])

    # No solutions were found.
    # TODO: Check if there exist some anyways?
    if len(partial_solutions) == 0:
        return []

    # Is this the last equation, that is, deepest hierarchy?
    if len(gb) == 1:
        return map(lambda s:[s], partial_solutions)

    # Finally, call this function recursively for each root replacing
    # the corresponding variable in the system.
    result = []
    for r in partial_solutions:
        new_system = []
        for eq in gb[:-1]:
            new_eq = eq.sympy_expr.subs(eq.var[-1], r).expand()
            if new_eq is not S.Zero:
                new_system.append(
                    Polynomial(new_eq, var=eq.var[:-1], order='lex'))
        if not new_system:
            return []
        for nps in equation_system(new_system):
            result.append(nps + [r])

    return result

def n_poly(f):
    """Checks if the polynomial can be simplifed by substituting the
    variable by a power.

    Returns a list of real or complex roots, or None, if no solution
    was found. Works for univariate instances of Polynomial only.
    """
    def roots_of_unity(n):
        result = []
        for i in range(0,n):
            # result.append((exp(2*i*pi*I/n)).evalc())
            result.append(exp(2*i*pi*I/n))
        return result

    exponents = map(lambda t:int(t[1]), f.coeffs)
    g = reduce(numbers.gcd, exponents)
    if g == 1:
        return None
    n = int(f.coeffs[0][1]/g)
    if not n in [1, 2, 3]: # Cases where solution can be computed
        return None

    ff = Polynomial(coeffs=tuple(map(lambda t:(t[0], t[1]/g), f.coeffs)),
                    var=f.var, order=f.order)
    # TODO: use uv() here?
    if n == 1:
        sol = [-(ff.coeffs[1][0]/ff.coeffs[0][0])]
    if n == 2:
        sol = quadratic(ff)
    elif n == 3:
        sol = cubic(ff)
    return [(zeta*s**Rational(1,g)).expand()
            for s in sol for zeta in roots_of_unity(g)]
        
def quadratic(f):
    """Returns the real or complex roots of a quadratic polynomial.

    Works for univariate instances of Polynomial only."""

    # Get monic polynomial, for p-q formula
    f = f.as_monic()[1]

    # Solve special cases:
    if len(f.coeffs) == 1:
        return [S.Zero]
    if len(f.coeffs) == 2:
        if f.coeffs[1][1] == 1: # No constant term
            return [S.Zero, -(f.coeffs[1][0])]
        else: # No linear term
            q = -(f.coeffs[1][0])
            if q > 0:
                return [-sqrt(q), sqrt(q)]
            else:
                return [-sqrt(q), sqrt(q)]
                
    p = f.coeffs[1][0]
    q = f.coeffs[2][0]
    discr = p**2 - 4*q
    if (not discr.is_complex) or discr > 0:
        return [-p/2 + sqrt(discr)/2,
                -p/2 - sqrt(discr)/2]
    elif discr == 0:
        return [-p/2]
    else: # discr < 0
        return [-p/2 + I*sqrt(-discr)/2,
                -p/2 - I*sqrt(-discr)/2]

# TODO: Implement function to find roots of quartic polynomials?

def rat_roots(f):
    """Returns a list of rational roots of an integer Polynomial.

    For an polynomial an*x**n + ... + a0, all rational roots are of
    the form p/q, where p and q are integer factors of a0 and an.
    """
    an_divs = integer_divisors(int(f.coeffs[0][0]))
    a0_divs = integer_divisors(int(f.coeffs[-1][0]))
    result = []
    for p in a0_divs:
        for q in an_divs:
            if f(Rational(p, q)) is S.Zero:
                result.append(Rational(p, q))
            if f(Rational(-p, q)) is S.Zero:
                result.append(Rational(-p, q))
    # Now check, if 0 is a root
    if f.sympy_expr.subs(f.var[0], S.Zero).expand() is S.Zero:
        result.append(S.Zero)
    return result

def real_roots(s, a=None, b=None):
    """Returns the number of unique real roots of f in the interval (a, b].

    Assumes a sturm sequence of an univariate, square-free instance of
    Polynomial with real coeffs. The boundaries a and b can be omitted
    to check the whole real line.
    """
    def sign_changes(list):
        counter = 0
        current = list[0]
        for el in list:
            if (current < 0 and el >= 0) or \
               (current > 0 and el <= 0):
                counter += 1
            current = el
        return counter
    
    if a is None: # a = -oo
        sa = sign_changes(map(
            lambda p:p.coeffs[0][0]*(-1)**p.coeffs[0][1], s))
    else:
        sa = sign_changes(map(lambda p:p.sympy_expr.subs(p.var[0], a), s))
    if b is None: # b = oo
        sb = sign_changes(map(lambda p:p.coeffs[0][0], s))
    else:
        sb = sign_changes(map(lambda p:p.sympy_expr.subs(p.var[0], b), s))
    return sa - sb
    
def sturm(f):
    """Compute the Sturm sequence of given polynomial

    Assumes an univariate instance of Polynomial with real coeffs.
    """
    seq = [f]
    seq.append(f.diff(f.var[0]))
    while seq[-1].sympy_expr is not S.Zero:
        seq.append(-(div_.mv(seq[-2], seq[-1])[-1]))
    return seq[:-1]

def uv(f, verbose=False):
    """Compute the roots of an univariate polynomial.
    """
    n = f.coeffs[0][1] # degree
    if n == 0: # constant
        return []
    elif n == 1:
        result = [-(f.coeffs[1][0] / f.coeffs[0][0])]
    elif n == 2:
        result = quadratic(f)
    elif n == 3:
        result = cubic(f)
    else:
        result = n_poly(f)
        if result is None:
            return []

    return result
