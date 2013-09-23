"""Algorithms for computing symbolic roots of polynomials. """

from __future__ import print_function, division

import math

from sympy.core.symbol import Dummy, Symbol, symbols
from sympy.core import S, I, pi
from sympy.core.mul import expand_2arg
from sympy.core.sympify import sympify
from sympy.core.numbers import Rational, igcd

from sympy.ntheory import divisors, isprime, nextprime
from sympy.functions import exp, sqrt, re, im, Abs, cos, sin

from sympy.polys.polytools import Poly, cancel, factor, gcd_list, discriminant
from sympy.polys.specialpolys import cyclotomic_poly
from sympy.polys.polyerrors import PolynomialError, GeneratorsNeeded, DomainError
from sympy.polys.polyquinticconst import PolyQuintic
from sympy.polys.rationaltools import together

from sympy.simplify import simplify, powsimp
from sympy.utilities import default_sort_key, public

from sympy.core.compatibility import reduce, xrange


def roots_linear(f):
    """Returns a list of roots of a linear polynomial."""
    r = -f.nth(0)/f.nth(1)
    dom = f.get_domain()

    if not dom.is_Numerical:
        if dom.is_Composite:
            r = factor(r)
        else:
            r = simplify(r)

    return [r]


def roots_quadratic(f):
    """Returns a list of roots of a quadratic polynomial."""
    a, b, c = f.all_coeffs()
    dom = f.get_domain()

    def _simplify(expr):
        if dom.is_Composite:
            return factor(expr)
        else:
            return simplify(expr)

    if c is S.Zero:
        r0, r1 = S.Zero, -b/a

        if not dom.is_Numerical:
            r1 = _simplify(r1)
    elif b is S.Zero:
        r = -c/a

        if not dom.is_Numerical:
            R = sqrt(_simplify(r))
        else:
            R = sqrt(r)

        r0 = R
        r1 = -R
    else:
        d = b**2 - 4*a*c

        if dom.is_Numerical:
            D = sqrt(d)

            r0 = (-b + D) / (2*a)
            r1 = (-b - D) / (2*a)
        else:
            D = sqrt(_simplify(d))
            A = 2*a

            E = _simplify(-b/A)
            F = D/A

            r0 = E + F
            r1 = E - F

    return sorted([expand_2arg(i) for i in (r0, r1)], key=default_sort_key)


def roots_cubic(f):
    """Returns a list of roots of a cubic polynomial."""
    _, a, b, c = f.monic().all_coeffs()

    if c is S.Zero:
        x1, x2 = roots([1, a, b], multiple=True)
        return [x1, S.Zero, x2]

    p = b - a**2/3
    q = c - a*b/3 + 2*a**3/27

    pon3 = p/3
    aon3 = a/3

    if p is S.Zero:
        if q is S.Zero:
            return [-aon3]*3
        else:
            if q.is_real:
                if (q > 0) is True:
                    u1 = -q**Rational(1, 3)
                else:
                    u1 = (-q)**Rational(1, 3)
            else:
                u1 = (-q)**Rational(1, 3)
    elif q is S.Zero:
        y1, y2 = roots([1, 0, p], multiple=True)
        return [tmp - aon3 for tmp in [y1, S.Zero, y2]]
    elif q.is_real and q < 0:
        u1 = -(-q/2 + sqrt(q**2/4 + pon3**3))**Rational(1, 3)
    else:
        u1 = (q/2 + sqrt(q**2/4 + pon3**3))**Rational(1, 3)

    coeff = S.ImaginaryUnit*sqrt(3)/2

    u2 = u1*(-S.Half + coeff)
    u3 = u1*(-S.Half - coeff)

    if p is S.Zero:
        return [u1 - aon3, u2 - aon3, u3 - aon3]

    soln = [
        -u1 + pon3/u1 - aon3,
        -u2 + pon3/u2 - aon3,
        -u3 + pon3/u3 - aon3
    ]

    return soln

def _roots_quartic_euler(p, q, r, a):
    """
    Descartes-Euler solution of the quartic equation

    Parameters
    ==========

    p, q, r: coefficients of ``x**4 + p*x**2 + q*x + r``
    a: shift of the roots

    Notes
    =====

    This is a helper function for ``roots_quartic``.

    Look for solutions of the form ::

      ``x1 = sqrt(R) - sqrt(A + B*sqrt(R))``
      ``x2 = -sqrt(R) - sqrt(A - B*sqrt(R))``
      ``x3 = -sqrt(R) + sqrt(A - B*sqrt(R))``
      ``x4 = sqrt(R) + sqrt(A + B*sqrt(R))``

    To satisfy the quartic equation one must have
    ``p = -2*(R + A); q = -4*B*R; r = (R - A)**2 - B**2*R``
    so that ``R`` must satisfy the Descartes-Euler resolvent equation
    ``64*R**3 + 32*p*R**2 + (4*p**2 - 16*r)*R - q**2 = 0``

    If the resolvent does not have a rational solution, return None;
    in that case it is likely that the Ferrari method gives a simpler
    solution.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.polys.polyroots import _roots_quartic_euler
    >>> p, q, r = -S(64)/5, -S(512)/125, -S(1024)/3125
    >>> _roots_quartic_euler(p, q, r, S(0))[0]
    -sqrt(32*sqrt(5)/125 + 16/5) + 4*sqrt(5)/5
    """
    from sympy.solvers import solve
    # solve the resolvent equation
    x = Symbol('x')
    eq = 64*x**3 + 32*p*x**2 + (4*p**2 - 16*r)*x - q**2
    xsols = list(roots(Poly(eq, x), cubics=False).keys())
    xsols = [sol for sol in xsols if sol.is_rational]
    if not xsols:
        return None
    R = max(xsols)
    c1 = sqrt(R)
    B = -q*c1/(4*R)
    A = -R - p/2
    c2 = sqrt(A + B)
    c3 = sqrt(A - B)
    return [c1 - c2 - a, -c1 - c3 - a, -c1 + c3 - a, c1 + c2 - a]


def roots_quartic(f):
    r"""
    Returns a list of roots of a quartic polynomial.

    There are many references for solving quartic expressions available [1-5].
    This reviewer has found that many of them require one to select from among
    2 or more possible sets of solutions and that some solutions work when one
    is searching for real roots but don't work when searching for complex roots
    (though this is not always stated clearly). The following routine has been
    tested and found to be correct for 0, 2 or 4 complex roots.

    The quasisymmetric case solution [6] looks for quartics that have the form
    `x**4 + A*x**3 + B*x**2 + C*x + D = 0` where `(C/A)**2 = D`.

    Although there is a general solution, simpler results can be obtained for
    certain values of the coefficients. In all cases, 4 roots are returned:

      1) `f = c + a*(a**2/8 - b/2) == 0`
      2) `g = d - a*(a*(3*a**2/256 - b/16) + c/4) = 0`
      3) if `f != 0` and `g != 0` and `p = -d + a*c/4 - b**2/12` then
        a) `p == 0`
        b) `p != 0`

    Examples
    ========

        >>> from sympy import Poly, symbols, I
        >>> from sympy.polys.polyroots import roots_quartic

        >>> r = roots_quartic(Poly('x**4-6*x**3+17*x**2-26*x+20'))

        >>> # 4 complex roots: 1+-I*sqrt(3), 2+-I
        >>> sorted(str(tmp.evalf(n=2)) for tmp in r)
        ['1.0 + 1.7*I', '1.0 - 1.7*I', '2.0 + 1.0*I', '2.0 - 1.0*I']

    References
    ==========

    1. http://mathforum.org/dr.math/faq/faq.cubic.equations.html
    2. http://en.wikipedia.org/wiki/Quartic_function#Summary_of_Ferrari.27s_method
    3. http://planetmath.org/encyclopedia/GaloisTheoreticDerivationOfTheQuarticFormula.html
    4. http://staff.bath.ac.uk/masjhd/JHD-CA.pdf
    5. http://www.albmath.org/files/Math_5713.pdf
    6. http://www.statemaster.com/encyclopedia/Quartic-equation
    7. eqworld.ipmnet.ru/en/solutions/ae/ae0108.pdf
    """
    _, a, b, c, d = f.monic().all_coeffs()

    if not d:
        return [S.Zero] + roots([1, a, b, c], multiple=True)
    elif (c/a)**2 == d:
        x, m = f.gen, c/a

        g = Poly(x**2 + a*x + b - 2*m, x)

        z1, z2 = roots_quadratic(g)

        h1 = Poly(x**2 - z1*x + m, x)
        h2 = Poly(x**2 - z2*x + m, x)

        r1 = roots_quadratic(h1)
        r2 = roots_quadratic(h2)

        return r1 + r2
    else:
        a2 = a**2
        e = b - 3*a2/8
        f = c + a*(a2/8 - b/2)
        g = d - a*(a*(3*a2/256 - b/16) + c/4)
        aon4 = a/4

        if f is S.Zero:
            y1, y2 = [sqrt(tmp) for tmp in
                      roots([1, e, g], multiple=True)]
            return [tmp - aon4 for tmp in [-y1, -y2, y1, y2]]
        if g is S.Zero:
            y = [S.Zero] + roots([1, 0, e, f], multiple=True)
            return [tmp - aon4 for tmp in y]
        else:
            # Descartes-Euler method, see [7]
            sols = _roots_quartic_euler(e, f, g, aon4)
            if sols:
                return sols
            # Ferrari method, see [1, 2]
            a2 = a**2
            e = b - 3*a2/8
            f = c + a*(a2/8 - b/2)
            g = d - a*(a*(3*a2/256 - b/16) + c/4)
            p = -e**2/12 - g
            q = -e**3/108 + e*g/3 - f**2/8
            TH = Rational(1, 3)
            if p.is_zero:
                y = -5*e/6 - q**TH
            elif p.is_number and p.is_comparable:
                # with p != 0 then u below is not 0
                root = sqrt(q**2/4 + p**3/27)
                r = -q/2 + root  # or -q/2 - root
                u = r**TH  # primary root of solve(x**3-r, x)
                y = -5*e/6 + u - p/u/3
            else:
                raise PolynomialError('cannot return general quartic solution')
            w = sqrt(e + 2*y)
            arg1 = 3*e + 2*y
            arg2 = 2*f/w
            ans = []
            for s in [-1, 1]:
                root = sqrt(-(arg1 + s*arg2))
                for t in [-1, 1]:
                    ans.append((s*w - t*root)/2 - aon4)

            return ans

def roots_binomial(f):
    """Returns a list of roots of a binomial polynomial."""
    n = f.degree()

    a, b = f.nth(n), f.nth(0)
    alpha = (-cancel(b/a))**Rational(1, n)

    if alpha.is_number:
        alpha = alpha.expand(complex=True)

    roots, I = [], S.ImaginaryUnit

    for k in xrange(n):
        zeta = exp(2*k*S.Pi*I/n).expand(complex=True)
        roots.append((alpha*zeta).expand(power_base=False))

    return sorted(roots, key=default_sort_key)


def _inv_totient_estimate(m):
    """
    Find ``(L, U)`` such that ``L <= phi^-1(m) <= U``.

    Examples
    ========

    >>> from sympy.polys.polyroots import _inv_totient_estimate

    >>> _inv_totient_estimate(192)
    (192, 840)
    >>> _inv_totient_estimate(400)
    (400, 1750)

    """
    primes = [ d + 1 for d in divisors(m) if isprime(d + 1) ]

    a, b = 1, 1

    for p in primes:
        a *= p
        b *= p - 1

    L = m
    U = int(math.ceil(m*(float(a)/b)))

    P = p = 2
    primes = []

    while P <= U:
        p = nextprime(p)
        primes.append(p)
        P *= p

    P //= p
    b = 1

    for p in primes[:-1]:
        b *= p - 1

    U = int(math.ceil(m*(float(P)/b)))

    return L, U


def roots_cyclotomic(f, factor=False):
    """Compute roots of cyclotomic polynomials. """
    L, U = _inv_totient_estimate(f.degree())

    for n in xrange(L, U + 1):
        g = cyclotomic_poly(n, f.gen, polys=True)

        if f == g:
            break
    else:  # pragma: no cover
        raise RuntimeError("failed to find index of a cyclotomic polynomial")

    roots = []

    if not factor:
        for k in xrange(1, n + 1):
            if igcd(k, n) == 1:
                roots.append(exp(2*k*S.Pi*I/n).expand(complex=True))
    else:
        g = Poly(f, extension=(-1)**Rational(1, n))

        for h, _ in g.factor_list()[1]:
            roots.append(-h.TC())

    return sorted(roots, key=default_sort_key)


def roots_quintic(f):
    """
    Calulate exact roots of a solvable quintic
    """
    result = []
    coeff_5, coeff_4, p, q, r, s = f.all_coeffs()

    # Eqn must be of the form x^5 + px^3 + qx^2 + rx + s
    if coeff_4:
        return result

    if coeff_5 != 1:
        l = [p/coeff_5, q/coeff_5, r/coeff_5, s/coeff_5]
        if not all(coeff.is_Rational for coeff in l):
            return result
        f = Poly(f/coeff_5)
    quintic = PolyQuintic(f)

    # Eqn standardised. Algo for solving starts here
    if not f.is_irreducible:
        return result

    f20 = quintic.f20
    # Check if f20 has linear factors over domain Z
    if f20.is_irreducible:
        return result

    # Now, we know that f is solvable
    for _factor in f20.factor_list()[1]:
        if _factor[0].is_linear:
            theta = _factor[0].root(0)
            break
    d = discriminant(f)
    delta = sqrt(d)
    # zeta = a fifth root of unity
    zeta1, zeta2, zeta3, zeta4 = quintic.zeta
    T = quintic.T(theta, d)
    tol = S(1e-10)
    alpha = T[1] + T[2]*delta
    alpha_bar = T[1] - T[2]*delta
    beta = T[3] + T[4]*delta
    beta_bar = T[3] - T[4]*delta

    disc = alpha**2 - 4*beta
    disc_bar = alpha_bar**2 - 4*beta_bar

    l0 = quintic.l0(theta)

    l1 = _quintic_simplify((-alpha + sqrt(disc)) / S(2))
    l4 = _quintic_simplify((-alpha - sqrt(disc)) / S(2))

    l2 = _quintic_simplify((-alpha_bar + sqrt(disc_bar)) / S(2))
    l3 = _quintic_simplify((-alpha_bar - sqrt(disc_bar)) / S(2))

    order = quintic.order(theta, d)
    test = (order*delta.n()) - ( (l1.n() - l4.n())*(l2.n() - l3.n()) )
    # Comparing floats
    # Problems importing on top
    from sympy.utilities.randtest import comp
    if not comp(test, 0, tol):
        l2, l3 = l3, l2

    # Now we have correct order of l's
    R1 = l0 + l1*zeta1 + l2*zeta2 + l3*zeta3 + l4*zeta4
    R2 = l0 + l3*zeta1 + l1*zeta2 + l4*zeta3 + l2*zeta4
    R3 = l0 + l2*zeta1 + l4*zeta2 + l1*zeta3 + l3*zeta4
    R4 = l0 + l4*zeta1 + l3*zeta2 + l2*zeta3 + l1*zeta4

    Res = [None, [None]*5, [None]*5, [None]*5, [None]*5]
    Res_n = [None, [None]*5, [None]*5, [None]*5, [None]*5]
    sol = Symbol('sol')

    # Simplifying improves performace a lot for exact expressions
    R1 = _quintic_simplify(R1)
    R2 = _quintic_simplify(R2)
    R3 = _quintic_simplify(R3)
    R4 = _quintic_simplify(R4)

    # Solve imported here. Causing problems if imported as 'solve'
    # and hence the changed name
    from sympy.solvers.solvers import solve as _solve
    a, b = symbols('a b', cls=Dummy)
    _sol = _solve( sol**5 - a - I*b, sol)
    for i in range(5):
        _sol[i] = factor(_sol[i])
    R1 = R1.as_real_imag()
    R2 = R2.as_real_imag()
    R3 = R3.as_real_imag()
    R4 = R4.as_real_imag()

    for i, root in enumerate(_sol):
        Res[1][i] = _quintic_simplify(root.subs({ a: R1[0], b: R1[1] }))
        Res[2][i] = _quintic_simplify(root.subs({ a: R2[0], b: R2[1] }))
        Res[3][i] = _quintic_simplify(root.subs({ a: R3[0], b: R3[1] }))
        Res[4][i] = _quintic_simplify(root.subs({ a: R4[0], b: R4[1] }))

    for i in range(1, 5):
        for j in range(5):
            Res_n[i][j] = Res[i][j].n()
            Res[i][j] = _quintic_simplify(Res[i][j])
    r1 = Res[1][0]
    r1_n = Res_n[1][0]

    for i in range(5):
        if comp(im(r1_n*Res_n[4][i]), 0, tol):
            r4 = Res[4][i]
            break

    u, v = quintic.uv(theta, d)
    sqrt5 = math.sqrt(5)

    # Now we have various Res values. Each will be a list of five
    # values. We have to pick one r value from those five for each Res
    u, v = quintic.uv(theta, d)
    testplus = (u + v*delta*sqrt(5)).n()
    testminus = (u - v*delta*sqrt(5)).n()

    # Evaluated numbers suffixed with _n
    # We will use evaluated numbers for calculation. Much faster.
    r4_n = r4.n()
    r2 = r3 = None

    for i in range(5):
        r2temp_n = Res_n[2][i]
        for j in range(5):
            # Again storing away the exact number and using
            # evaluated numbers in computations
            r3temp_n = Res_n[3][j]

            if( comp( r1_n*r2temp_n**2 + r4_n*r3temp_n**2 - testplus, 0, tol) and
                comp( r3temp_n*r1_n**2 + r2temp_n*r4_n**2 - testminus, 0, tol ) ):
                r2 = Res[2][i]
                r3 = Res[3][j]
                break
        if r2:
            break

    # Now, we have r's so we can get roots
    x1 = (r1 + r2 + r3 + r4)/5
    x2 = (r1*zeta4 + r2*zeta3 + r3*zeta2 + r4*zeta1)/5
    x3 = (r1*zeta3 + r2*zeta1 + r3*zeta4 + r4*zeta2)/5
    x4 = (r1*zeta2 + r2*zeta4 + r3*zeta1 + r4*zeta3)/5
    x5 = (r1*zeta1 + r2*zeta2 + r3*zeta3 + r4*zeta4)/5
    result = [x1, x2, x3, x4, x5]

    # Now check if solutions are distinct

    result_n = []
    for root in result:
        result_n.append(root.n(5))
    result_n = sorted(result_n, key=default_sort_key)

    prev_entry = None
    for r in result_n:
        if r == prev_entry:
            # Roots are identical. Abort. Return []
            # and fall back to usual solve
            return []
        prev_entry = r

    return result


def _quintic_simplify(expr):
    expr = powsimp(expr)
    expr = cancel(expr)
    return together(expr)


def _integer_basis(poly):
    """Compute coefficient basis for a polynomial over integers.

    Returns the integer ``div`` such that substituting ``x = div*y``
    ``p(x) = m*q(y)`` where the coefficients of ``q`` are smaller
    than those of ``p``.

    For example ``x**5 + 512*x + 1024 = 0``
    with ``div = 4`` becomes ``y**5 + 2*y + 1 = 0``

    Returns the integer ``div`` or ``None`` if there is no possible scaling.

    Examples
    ========

    >>> from sympy.polys import Poly
    >>> from sympy.abc import x
    >>> from sympy.polys.polyroots import _integer_basis
    >>> p = Poly(x**5 + 512*x + 1024, x, domain='ZZ')
    >>> _integer_basis(p)
    4
    """
    monoms, coeffs = list(zip(*poly.terms()))

    monoms, = list(zip(*monoms))
    coeffs = list(map(abs, coeffs))

    if coeffs[0] < coeffs[-1]:
        coeffs = list(reversed(coeffs))
        n = monoms[0]
        monoms = [n - i for i in reversed(monoms)]
    else:
        return None

    monoms = monoms[:-1]
    coeffs = coeffs[:-1]

    divs = reversed(divisors(gcd_list(coeffs))[1:])

    try:
        div = next(divs)
    except StopIteration:
        return None

    while True:
        for monom, coeff in zip(monoms, coeffs):
            if coeff % div**monom != 0:
                try:
                    div = next(divs)
                except StopIteration:
                    return None
                else:
                    break
        else:
            return div


def preprocess_roots(poly):
    """Try to get rid of symbolic coefficients from ``poly``. """
    coeff = S.One

    try:
        _, poly = poly.clear_denoms(convert=True)
    except DomainError:
        return coeff, poly

    poly = poly.primitive()[1]
    poly = poly.retract()

    # TODO: This is fragile. Figure out how to make this independent of construct_domain().
    if poly.get_domain().is_Poly and all(c.is_term for c in poly.rep.coeffs()):
        poly = poly.inject()

        strips = list(zip(*poly.monoms()))
        gens = list(poly.gens[1:])

        base, strips = strips[0], strips[1:]

        for gen, strip in zip(list(gens), strips):
            reverse = False

            if strip[0] < strip[-1]:
                strip = reversed(strip)
                reverse = True

            ratio = None

            for a, b in zip(base, strip):
                if not a and not b:
                    continue
                elif not a or not b:
                    break
                elif b % a != 0:
                    break
                else:
                    _ratio = b // a

                    if ratio is None:
                        ratio = _ratio
                    elif ratio != _ratio:
                        break
            else:
                if reverse:
                    ratio = -ratio

                poly = poly.eval(gen, 1)
                coeff *= gen**(-ratio)
                gens.remove(gen)

        if gens:
            poly = poly.eject(*gens)

    if poly.is_univariate and poly.get_domain().is_ZZ:
        basis = _integer_basis(poly)

        if basis is not None:
            n = poly.degree()

            def func(k, coeff):
                return coeff//basis**(n - k[0])

            poly = poly.termwise(func)
            coeff *= basis

    return coeff, poly


@public
def roots(f, *gens, **flags):
    """
    Computes symbolic roots of a univariate polynomial.

    Given a univariate polynomial f with symbolic coefficients (or
    a list of the polynomial's coefficients), returns a dictionary
    with its roots and their multiplicities.

    Only roots expressible via radicals will be returned.  To get
    a complete set of roots use RootOf class or numerical methods
    instead. By default cubic and quartic formulas are used in
    the algorithm. To disable them because of unreadable output
    set ``cubics=False`` or ``quartics=False`` respectively.

    To get roots from a specific domain set the ``filter`` flag with
    one of the following specifiers: Z, Q, R, I, C. By default all
    roots are returned (this is equivalent to setting ``filter='C'``).

    By default a dictionary is returned giving a compact result in
    case of multiple roots.  However to get a tuple containing all
    those roots set the ``multiple`` flag to True.

    Examples
    ========

    >>> from sympy import Poly, roots
    >>> from sympy.abc import x, y

    >>> roots(x**2 - 1, x)
    {-1: 1, 1: 1}

    >>> p = Poly(x**2-1, x)
    >>> roots(p)
    {-1: 1, 1: 1}

    >>> p = Poly(x**2-y, x, y)

    >>> roots(Poly(p, x))
    {-sqrt(y): 1, sqrt(y): 1}

    >>> roots(x**2 - y, x)
    {-sqrt(y): 1, sqrt(y): 1}

    >>> roots([1, 0, -1])
    {-1: 1, 1: 1}

    """
    from sympy.polys.polytools import to_rational_coeffs
    flags = dict(flags)

    auto = flags.pop('auto', True)
    cubics = flags.pop('cubics', True)
    quartics = flags.pop('quartics', True)
    quintics = flags.pop('quintics', False)
    multiple = flags.pop('multiple', False)
    filter = flags.pop('filter', None)
    predicate = flags.pop('predicate', None)

    if isinstance(f, list):
        if gens:
            raise ValueError('redundant generators given')

        x = Dummy('x')

        poly, i = {}, len(f) - 1

        for coeff in f:
            poly[i], i = sympify(coeff), i - 1

        f = Poly(poly, x, field=True)
    else:
        try:
            f = Poly(f, *gens, **flags)
        except GeneratorsNeeded:
            if multiple:
                return []
            else:
                return {}

        if f.is_multivariate:
            raise PolynomialError('multivariate polynomials are not supported')

    def _update_dict(result, root, k):
        if root in result:
            result[root] += k
        else:
            result[root] = k

    def _try_decompose(f):
        """Find roots using functional decomposition. """
        factors, roots = f.decompose(), []

        for root in _try_heuristics(factors[0]):
            roots.append(root)

        for factor in factors[1:]:
            previous, roots = list(roots), []

            for root in previous:
                g = factor - Poly(root, f.gen)

                for root in _try_heuristics(g):
                    roots.append(root)

        return roots

    def _try_heuristics(f):
        """Find roots using formulas and some tricks. """
        if f.is_ground:
            return []
        if f.is_monomial:
            return [S(0)]*f.degree()

        if f.length() == 2:
            if f.degree() == 1:
                return list(map(cancel, roots_linear(f)))
            else:
                return roots_binomial(f)

        result = []

        for i in [-1, 1]:
            if not f.eval(i):
                f = f.quo(Poly(f.gen - i, f.gen))
                result.append(i)
                break

        n = f.degree()

        if n == 1:
            result += list(map(cancel, roots_linear(f)))
        elif n == 2:
            result += list(map(cancel, roots_quadratic(f)))
        elif f.is_cyclotomic:
            result += roots_cyclotomic(f)
        elif n == 3 and cubics:
            result += roots_cubic(f)
        elif n == 4 and quartics:
            result += roots_quartic(f)
        elif n == 5 and quintics:
            result += roots_quintic(f)

        return result

    (k,), f = f.terms_gcd()

    if not k:
        zeros = {}
    else:
        zeros = {S(0): k}

    coeff, f = preprocess_roots(f)

    if auto and f.get_domain().has_Ring:
        f = f.to_field()

    rescale_x = None
    translate_x = None

    result = {}

    if not f.is_ground:
        if not f.get_domain().is_Exact:
            for r in f.nroots():
                _update_dict(result, r, 1)
        elif f.degree() == 1:
            result[roots_linear(f)[0]] = 1
        elif f.degree() == 2:
            for r in roots_quadratic(f):
                _update_dict(result, r, 1)
        elif f.length() == 2:
            for r in roots_binomial(f):
                _update_dict(result, r, 1)
        else:
            _, factors = Poly(f.as_expr()).factor_list()

            if len(factors) == 1 and factors[0][1] == 1:
                if f.get_domain().is_EX:
                    res = to_rational_coeffs(f)
                    if res:
                        if res[0] is None:
                            translate_x, f = res[2:]
                        else:
                            rescale_x, f = res[1], res[-1]
                        result = roots(f)
                        if not result:
                            for root in _try_decompose(f):
                                _update_dict(result, root, 1)
                else:
                    for root in _try_decompose(f):
                        _update_dict(result, root, 1)
            else:
                for factor, k in factors:
                    for r in _try_heuristics(Poly(factor, f.gen, field=True)):
                        _update_dict(result, r, k)

    if coeff is not S.One:
        _result, result, = result, {}

        for root, k in _result.items():
            result[coeff*root] = k

    result.update(zeros)

    if filter not in [None, 'C']:
        handlers = {
            'Z': lambda r: r.is_Integer,
            'Q': lambda r: r.is_Rational,
            'R': lambda r: r.is_real,
            'I': lambda r: r.is_imaginary,
        }

        try:
            query = handlers[filter]
        except KeyError:
            raise ValueError("Invalid filter: %s" % filter)

        for zero in dict(result).keys():
            if not query(zero):
                del result[zero]

    if predicate is not None:
        for zero in dict(result).keys():
            if not predicate(zero):
                del result[zero]
    if rescale_x:
        result1 = {}
        for k, v in result.items():
            result1[k*rescale_x] = v
        result = result1
    if translate_x:
        result1 = {}
        for k, v in result.items():
            result1[k + translate_x] = v
        result = result1

    if not multiple:
        return result
    else:
        zeros = []

        for zero, k in result.items():
            zeros.extend([zero]*k)

        return sorted(zeros, key=default_sort_key)


def root_factors(f, *gens, **args):
    """
    Returns all factors of a univariate polynomial.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.polys.polyroots import root_factors

    >>> root_factors(x**2 - y, x)
    [x - sqrt(y), x + sqrt(y)]

    """
    args = dict(args)
    filter = args.pop('filter', None)

    F = Poly(f, *gens, **args)

    if not F.is_Poly:
        return [f]

    if F.is_multivariate:
        raise ValueError('multivariate polynomials not supported')

    x = F.gens[0]

    zeros = roots(F, filter=filter)

    if not zeros:
        factors = [F]
    else:
        factors, N = [], 0

        for r, n in zeros.items():
            factors, N = factors + [Poly(x - r, x)]*n, N + n

        if N < F.degree():
            G = reduce(lambda p, q: p*q, factors)
            factors.append(F.quo(G))

    if not isinstance(f, Poly):
        factors = [ f.as_expr() for f in factors ]

    return sorted(factors, key=default_sort_key)
