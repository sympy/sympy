
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.basic import Basic, S
from sympy.core.sympify import sympify
from sympy.core.numbers import Rational
from sympy.core.methods import NoRelMeths, ArithMeths

from sympy.polys.polynomial import Poly, PolynomialError
from sympy.polys.algorithms import poly_decompose, poly_sqf, poly_div

from sympy.functions import exp

def roots_linear(f):
    """Returns a list of roots of a linear polynomial."""
    return [-f.coeff(0)/f.coeff(1)]

def roots_quadratic(f):
    """Returns a list of roots of a quadratic polynomial."""
    a, b, c = f.iter_all_coeffs()

    if c is S.Zero:
        return [c, -b/a]

    d = (b**2 - 4*a*c)**S.Half

    roots = [
        (-b + d) / (2*a),
        (-b - d) / (2*a),
    ]

    return [ r.expand() for r in roots ]

def roots_cubic(f):
    """Returns a list of  roots of a cubic polynomial."""
    a, b, c, d = f.iter_all_coeffs()

    b, c, d = b/a, c/a, d/a

    p = c - b**2 / 3
    q = d + (2*b**3 - 9*b*c) / 27

    if p is S.Zero:
        if q is S.Zero:
            return [-b/3] * 3
        else:
            u1 = q**Rational(1, 3)
    else:
        u1 = (q/2 + (q**2/4 + p**3/27)**S.Half)**Rational(1, 3)

    coeff = S.ImaginaryUnit*3**S.Half / 2

    u2 = u1*(-S.Half + coeff)
    u3 = u1*(-S.Half - coeff)

    roots = [
        (p/(3*u1) - u1 - b/3),
        (p/(3*u2) - u2 - b/3),
        (p/(3*u3) - u3 - b/3),
    ]

    return [ r.expand() for r in roots ]

def roots_quartic(f):
    """Returns a list of roots of a quartic polynomial."""
    a, b, c, d, e = f.iter_all_coeffs()

    p = -2*b
    q = b**2 + a*c - 4*d
    r = c**2 + a**2*d - a*b*c

    r = cubic(1, p, q, r)

    u = r[1] + r[2] - r[0]
    v = (u**2 - 16*d)**S.Half
    w = (a**2 - 4*r[0])**S.Half

    A = (u + v) / 4
    B = (u - v) / 4
    C = (-a + w) / 2
    D = (-a - w) / 2

    E = (C**2 - 4*A)**S.Half
    F = (D**2 - 4*B)**S.Half

    roots = [
        (C + E) / 2,
        (C - E) / 2,
        (D + F) / 2,
        (D - F) / 2,
    ]

    return [ r.expand() for r in roots ]

def roots_binomial(f):
    """Returns a list of roots of a binomial polynomial."""
    n = f.degree

    a, b = f.coeff(n), f.coeff(0)
    alpha = (-b/a)**Rational(1, n)

    if alpha.is_number:
        alpha = alpha.expand(complex=True)

    roots, I = [], S.ImaginaryUnit

    for k in xrange(n):
        zeta = exp(2*k*S.Pi*I/n).expand(complex=True)
        roots.append((alpha*zeta).expand())

    return roots

def roots(f, *symbols, **flags):
    """Computes symbolic roots of an univariate polynomial.

       Given an univariate  polynomial f with symbolic coefficients,
       returns a dictionary with its roots and their multiplicities.

       Only roots expressible via radicals will be returned.  To get
       a complete set of roots use RootOf class or numerical methods
       instead. By default cubic and quartic formulas aren't used in
       the algorithm because of unreadable output. To enable them
       set cubics=True or quartics=True respectively.

       To get roots from a specific domain set the 'domain' flag with
       one of the following specifiers: Z, Q, R, I, C. By default all
       roots are returned (this is equivalent to setting domain='C').

       By default a dictionary is returned giving a compact result in
       case of multiple roots.  However to get a tuple containing all
       those roots set the 'multiple' flag to True.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> roots(x**2 - 1, x)
       {1: 1, -1: 1}

       >>> roots(x**2 - y, x)
       {-y**(1/2): 1, y**(1/2): 1}

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise PolynomialError

    if f.is_multivariate:
        raise PolynomialError

    def roots_trivial(g):
        if g.length == 1:
            if g.is_constant:
                return []
            else:
                return [S.Zero] * g.degree

        (k,), g = g.as_reduced()

        if k == 0:
            zeros = []
        else:
            zeros = [S.Zero] * k

        if g.length == 2:
            zeros += roots_binomial(g)
        else:
            n = g.degree

            if n == 2:
                zeros += roots_quadratic(g)
            else:
                if n == 3 and flags.get('cubics', False):
                    from sympy.polynomials.factor_ import factor
                    from sympy.polynomials.base import PolynomialException

                    try:
                        factors = factor(g.as_basic(), g.symbols)

                        if len(factors) <= 2:
                            raise PolynomialException

                        for term in factors:
                            zeros += roots(term.as_poly(*g.symbols))
                    except PolynomialException:
                        zeros += roots_cubic(g)
                elif n == 4 and flags.get('quartics', False):
                    zeros += roots_quartic(g)

        return zeros

    multiple = flags.get('multiple', False)

    if f.length == 1:
        if f.is_constant:
            if multiple:
                return ()
            else:
                return {}
        else:
            if multiple:
                return (S.Zero,) * f.degree
            else:
                return { S.Zero : f.degree }

    (k,), f = f.as_reduced()

    if k == 0:
        result = {}
    else:
        result = { S.Zero : k }

    factors = poly_decompose(f)

    zeros, g = {}, factors[0]

    for i, h in enumerate(poly_sqf(g)):
        for zero in roots_trivial(h):
            if zeros.has_key(zero):
                zeros[zero] += i+1
            else:
                zeros[zero] = i+1

    for factor in factors[1:]:
        previous, zeros = zeros.copy(), {}

        for zero, i in previous.iteritems():
            g = factor.sub_term(zero, (0,))

            for j, h in enumerate(poly_sqf(g)):
                for zero in roots_trivial(h):
                    if zeros.has_key(zero):
                        zeros[zero] += i*(j+1)
                    else:
                        zeros[zero] = i*(j+1)

    domain = flags.get('domain', None)

    if domain not in [None, 'C']:
        handlers = {
            'Z' : lambda r: r.is_Integer,
            'Q' : lambda r: r.is_Rational,
            'R' : lambda r: r.is_real,
            'I' : lambda r: r.is_imaginary,
        }

        try:
            query = handlers[domain]
        except KeyError:
            raise ValueError("Invalid domain: " + domain)

        for zero in dict(zeros).iterkeys():
            if not query(zero):
                del zeros[zero]

    result.update(zeros)

    if not multiple:
        return result
    else:
        zeros = ()

        for zero, k in result.iteritems():
            zeros += (zero,) * k

        return zeros

def poly_factors(f, *symbols, **flags):
    """Returns all factors of an univariate polynomial.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> factors = poly_factors(x**2-y, x)

       >>> factors[0].as_basic()
       x + y**(1/2)
       >>> factors[1].as_basic()
       x - y**(1/2)

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise PolynomialError

    if f.is_multivariate:
        raise PolynomialError
    else:
        x = f.symbols[0]

    zeros = roots(f, **flags)

    if not zeros:
        return (f,)
    else:
        factors, N = [], 0

        for r, n in zeros.iteritems():
            h = Poly([(S.One, 1), (-r, 0)], x)
            factors, N = factors + [h]*n, N + n

        if N < f.degree:
            g = reduce(lambda p,q: p*q, factors)
            factors.append(poly_div(f, g)[0])

        return tuple(factors)

def poly_sturm(f, *symbols):
    """Computes the Sturm sequence of a given polynomial.

       Given an univariate, square-free polynomial f(x) returns an
       associated Sturm sequence f_0(x), ..., f_n(x) defined by:

           f_0(x), f_1(x) = f(x), f'(x)
           f_n = -rem(f_{n-2}(x), f_{n-1}(x))

       For more information on the implemented algorithm refer to:

       [1] J.H. Davenport, Y. Siret, E. Tournier, Computer Algebra
           Systems and Algorithms for Algebraic Computation,
           Academic Press, London, 1988, pp. 124-128

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise PolynomialError

    if f.is_multivariate:
        raise PolynomialError
    else:
        f = f.as_squarefree()

    sturm = [f, f.diff(f.symbols[0])]

    while not sturm[-1].is_zero:
        sturm.append(-poly_div(sturm[-2], sturm[-1])[1])

    return tuple(sturm[:-1])

_exact_roots_cache = {}

def _exact_roots(f):
    if _exact_roots_cache.has_key(f):
        zeros = _exact_roots_cache[f]
    else:
        exact, zeros = roots(f), []

        for zero, k in exact.iteritems():
            zeros += [zero] * k

        _exact_roots_cache[f] = zeros

    return zeros

class RootOf(Basic, NoRelMeths, ArithMeths):
    """Represents n-th root of an univariate polynomial. """

    def __new__(cls, f, index):
        if isinstance(f, RootsOf):
            f = f.poly
        elif not isinstance(f, Poly):
            raise PolynomialError

        if f.is_multivariate:
            raise PolynomialError

        if index < 0 or index >= f.degree:
            raise IndexError, "Index must be in [0, %d] range" % (f.degree-1)
        else:
            exact = _exact_roots(f)

            if index < len(exact):
                return exact[index]
            else:
                return Basic.__new__(cls, f, index)

    @property
    def poly(self):
        return self._args[0]

    @property
    def index(self):
        return self._args[1]

    def atoms(self, *args, **kwargs):
        return self.poly.atoms(*args, **kwargs)

class RootsOf(Basic, NoRelMeths, ArithMeths):
    """Represents all roots of an univariate polynomial.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> roots = RootsOf(x**2 + 2, x)

       >>> list(roots.roots())
       [I*2**(1/2), -I*2**(1/2)]

    """

    def __new__(cls, f, x=None):
        if not isinstance(f, Poly):
            f = Poly(f, x)
        elif x is not None:
            raise PolynomialError

        if f.is_multivariate:
            raise PolynomialError

        return Basic.__new__(cls, f)

    @property
    def poly(self):
        return self._args[0]

    @property
    def count(self):
        return self.poly.degree

    def roots(self):
        exact = _exact_roots(self.poly)

        for root in exact:
            yield root

        for j in range(len(exact), self.count):
            yield RootOf(self.poly, j)

    def exact_roots(self):
        exact = _exact_roots(self.poly)

        for root in exact:
            yield root

    def formal_roots(self):
        exact = _exact_roots(self.poly)

        for j in range(len(exact), self.count):
            yield RootOf(self.poly, j)

    def __call__(self, index):
        return RootOf(self.poly, index)

    def atoms(self, *args, **kwargs):
        return self.poly.atoms(*args, **kwargs)

class RootSum(Basic, NoRelMeths, ArithMeths):
    """Represents a sum of all roots of an univariate polynomial. """

    def __new__(cls, f, *args):
        roots = RootsOf(*args)

        if roots.count == 0:
            return S.Zero
        else:
            result = []

            for root in roots.exact_roots():
                result.append(f(root))

            if len(result) < roots.count:
                result.append(Basic.__new__(cls, f, roots))

            return Add(*result)

    @property
    def function(self):
        return self._args[0]

    @property
    def roots(self):
        return self._args[1]

    def doit(self, **hints):
        if not hints.get('roots', True):
            return self
        else:
            result = S.Zero

            for root in self.roots.formal_roots():
                result += self.function(root)

            return result
