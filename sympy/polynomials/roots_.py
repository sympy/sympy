"""Algorithms to determine the roots of polynomials"""

from sympy.polynomials.base import *
from sympy.polynomials import div_, groebner_

def cubic(f):
    """Computes the roots of a cubic polynomial.

    Usage:
    ======
        This function is called by the wrapper L{roots}, don't use it
        directly. The input is assumed to be a univariate instance of
        Polynomial of degree 3.

    References:
    ===========
        http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method

    """

    # Get monic polynomial.
    f = f.as_monic()[1]

    a = f.nth_coeff(2)
    b = f.nth_coeff(1)
    c = f.nth_coeff(0)

    # Substitute variable to get depressed cubic: t**3 + p * t + q
    p = b + -a**2/3
    q = c + (2*a**3 - 9*a*b)/27

    if p is S.Zero: # Solve special cases:
        if q is S.Zero:
            return [-a/3]
        else:
            u1 = q**Rational(1, 3)
    else:
        u1 = (q/2 + S.Sqrt(q**2/4 + p**3/27))**Rational(1, 3)

    u2 = u1*(Rational(-1, 2) + S.ImaginaryUnit*S.Sqrt(3)/2)
    u3 = u1*(Rational(-1, 2) - S.ImaginaryUnit*S.Sqrt(3)/2)

    return map(lambda u: (p/(u*3) - u - a/3).expand(), [u1, u2, u3])


def n_poly(f):
    """Tries to substitute a power of a variable, to simplify.

    Usage:
    ======
        This function is called by the wrapper L{roots}, don't use it
        directly. It returns 'None' if no such simplifcation is
        possible. The input f is assumed to be a univariate instance
        of Polynomial.

    References:
    ===========
        http://en.wikipedia.org/wiki/Root_of_unity
        http://en.wikipedia.org/wiki/Radical_root#Positive_real_numbers

    """

    def roots_of_unity(n):
        """Computes the list of the n-th roots of unity."""
        result = []
        for i in range(0,n):
            result.append(S.Exp(2*i*S.Pi*S.ImaginaryUnit/n))
        return result

    exponents = map(lambda t:int(t[1]), f.coeffs)
    g = reduce(numbers.gcd, exponents)
    if g == 1 or g == 0:
        return None
    n = int(f.coeffs[0][1]/g)
    if not n in [1, 2, 3]: # Cases where solution can be computed
        return None

    ff = Polynomial(coeffs=tuple(map(lambda t:(t[0], t[1]/g), f.coeffs)),
                    var=f.var, order=f.order)

    return [(zeta*s**Rational(1,g)).expand(complex=True)
            for s in roots(ff) for zeta in roots_of_unity(g)]


def quadratic(f):
    """Computes the roots of a quadratic polynomial.

    Usage:
    ======
        This function is called by the wrapper L{roots}, don't use it
        directly. The input is assumed to be a univariate instance of
        Polynomial of degree 2.

    References:
    ===========
        http://en.wikipedia.org/wiki/Quadratic_equation#Quadratic_formula

    """

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
                return [-S.Sqrt(q), S.Sqrt(q)]
            else:
                return [-S.Sqrt(q), S.Sqrt(q)]

    p = f.coeffs[1][0]
    q = f.coeffs[2][0]
    discr = p**2 - 4*q
    if (not discr.is_complex) or discr > 0:
        return [-p/2 + S.Sqrt(discr)/2,
                -p/2 - S.Sqrt(discr)/2]
    elif discr == 0:
        return [-p/2]
    else: # discr < 0
        return [-p/2 + S.ImaginaryUnit*S.Sqrt(-discr)/2,
                -p/2 - S.ImaginaryUnit*S.Sqrt(-discr)/2]


# TODO: Implement function to find roots of quartic polynomials?


def rat_roots(f):
    """Computes the rational roots of a polynomial.

    Usage:
    ======
        This function is called by the wrapper L{roots}, don't use it
        directly. The input is assumed to be a univariate and
        square-free instance of Polynomial, with integer coefficients.

    References:
    ===========
        http://en.wikipedia.org/wiki/Rational_root_theorem

    """

    # For an integer polynomial an*x**n + ... + a0, all rational roots
    # are of the form p/q, where p and q are factors of a0 and an.
    an_divs = integer_divisors(int(f.coeffs[0][0]))
    a0_divs = integer_divisors(int(f.coeffs[-1][0]))
    result = []
    for p in a0_divs:
        for q in an_divs:
            if f(Rational(p, q)) is S.Zero:
                result.append(Rational(p, q))
            if f(Rational(-p, q)) is S.Zero:
                result.append(Rational(-p, q))
    # Now check if 0 is a root.
    if f.sympy_expr.subs(f.var[0], S.Zero).expand() is S.Zero:
        result.append(S.Zero)
    return result

def count_real_roots(s, a=None, b=None):
    """Returns the number of unique real roots of f in the interval (a, b].

    Usage:
    ======
        The input can be a square-free and univariate polynomial, or a
        precomputed Sturm sequence, if you want to check one specific
        polynomial with several intervals. See L{sturm}.

        The boundaries a and b can be omitted to check the whole real
        line or one ray.

    Examples:
    =========
        >>> x = Symbol('x')
        >>> count_real_roots(x**2 - 1)
        2
        >>> count_real_roots(x**2 - 1, 0, 2)
        1

    References:
    ===========
        Davenport, Siret, Tournier: Computer Algebra, 1988

    """

    def sign_changes(lisp):
        """Counts how often the sign of consecutive list elements"""
        counter = 0
        current = lisp[0]
        for el in lisp:
            if (current < 0 and el >= 0) or \
               (current > 0 and el <= 0):
                counter += 1
            current = el
        return counter

    # Allow a polynomial instead of its Sturm sequence
    if not isinstance(s, list):
        s = sturm(s)

    if a is not None:
        a = sympify(a)
    if b is not None:
        b = sympify(b)

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
    """Compute the Sturm sequence of given polynomial.

    Usage:
    ======
        The input is assumed to be a square-free and univariate
        polynomial, either as a SymPy expression or as instance of
        Polynomial.

        The output is a list representing f's Sturm sequence, which is
        built similarly to the euclidian algorithm, beginning with f
        and its derivative.

        The result can be used in L{count_real_roots}.

    References:
    ===========
        Davenport, Siret, Tournier: Computer Algebra, 1988

    """

    if not isinstance(f, Polynomial):
        f = Polynomial(f)
    seq = [f]
    seq.append(f.diff(f.var[0]))
    while seq[-1].sympy_expr is not S.Zero:
        seq.append(-(div_.div(seq[-2], seq[-1])[-1]))
    return seq[:-1]


def roots(f, var=None):
    """Compute the roots of a univariate polynomial.

    Usage:
    ======
        The input f is assumed to be a univariate polynomial, either
        as SymPy expression or as instance of Polynomial. In the
        latter case, you can optionally specify the variable with
        'var'.

        The output is a list of all found roots with multiplicity.

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> roots(x**2 - 1)
        [1, -1]
        >>> roots(x - y, x)
        [y]

    Also see L{factor_.factor}, L{quadratic}, L{cubic}. L{n-poly},
    L{count_real_roots}.

    """

    from sympy.polynomials import factor_

    if not isinstance(f, Polynomial):
        f = Polynomial(f, var=var, order=None)
    if len(f.var) == 0:
        return []
    if len(f.var) > 1:
        raise PolynomialException('Multivariate polynomials not supported.')

    # Determine type of coeffs (for factorization purposes)
    symbols = f.sympy_expr.atoms(type=Symbol)
    symbols = filter(lambda a: not a in f.var, symbols)
    if symbols:
        coeff = 'sym'
    else:
        coeff = coeff_ring(get_numbers(f.sympy_expr))

    if coeff == 'rat':
        denom, f = f.as_integer()
        coeff = 'int'
    if coeff == 'int':
        content, f = f.as_primitive()
        # Hack to get some additional cases right:
        result = n_poly(f)
        if result is not None:
            return result
        factors = factor_.factor(f)
    else: # It's not possible to factorize.
        factors = [f]

    # Now check for roots in each factor.
    result = []
    for p in factors:
        n = p.coeffs[0][1] # Degree of the factor.
        if n == 0: # We have a constant factor.
            pass
        elif n == 1:
            if len(p.coeffs) == 2:
                result += [-(p.coeffs[1][0] / p.coeffs[0][0])]
            else:
                result += [S.Zero]
        elif n == 2:
            result += quadratic(p)
        elif n == 3:
            result += cubic(p)
        else:
            res = n_poly(p)
            if res is not None:
                result += res

    # With symbols, __nonzero__ returns a StrictInequality, Exception.
    try: result.sort()
    except: pass

    return result


def solve_system(eqs, var=None, order=None):
    """Solves a system of polynomial equations.

    Usage:
    ======
        Assumes to get a list of polynomials, either as SymPy
        expressions or instances of Polynomial. In the first case, you
        should specify the variables and monomial order through 'var'
        and 'order'. Otherwise, the polynomials should have matching
        variables and orders already. Only the first polynomial is
        checked for its type.

        This algorithm uses variable elimination and only works for
        zero-dimensional varieties, that is, a finite number of
        solutions, which is currently not tested.

    Examples:
    =========
        >>> x, y = symbols('xy')
        >>> f = y - x
        >>> g = x**2 + y**2 - 1
        >>> solve_system([f, g])
        [(-(1/2)**(1/2), -(1/2)**(1/2)), ((1/2)**(1/2), (1/2)**(1/2))]

    References:
    ===========
        Cox, Little, O'Shea: Ideals, Varieties and Algorithms,
        Springer, 2. edition, p. 113

    """

    def is_uv(f):
        """Is an instance of Polynomial univariate in its last variable?"""
        for term in f.coeffs:
            for exponent in term[1:-1]:
                if exponent > 0:
                    return False
        return True

    if not isinstance(eqs, list):
        eqs = [eqs]
    if not isinstance(eqs[0], Polynomial):
        if var is None:
            var = merge_var(*[f.atoms(type=Symbol) for f in eqs])
        eqs = [Polynomial(f, var=var, order='lex') for f in eqs]
    else:
        eqs = [Polynomial(f.sympy_expr, var=f.var, order='lex') for f in eqs]

    # First compute a Groebner base with the polynomials,
    # with lexicographic ordering, so that the last polynomial is
    # univariate and can be solved.
    gb = groebner_.groebner(eqs)

    # Now filter the the base elements, to get only the univariate ones.
    eliminated = filter(is_uv, gb)
    if len(eliminated) != 1:
        raise PolynomialException("System currently not solvable.")

    # Try to solve the polynomials with var eliminated.
    f = eliminated[0]
    partial_solutions = set(roots(f.sympy_expr, var=f.var[-1]))

    # No solutions were found.
    # TODO: Check if there exist some anyways?
    if len(partial_solutions) == 0:
        return []

    # Is this the last equation, that is, deepest hierarchy?
    if len(gb) == 1:
        return map(lambda s:(s,), partial_solutions)

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
        for nps in solve_system(new_system):
            result.append(nps + (r,))

    # With symbols, __nonzero__ returns a StrictInequality, Exception.
    try: result.sort()
    except: pass

    return result
