
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.basic import Basic, S
from sympy.core.sympify import sympify
from sympy.core.numbers import Rational

from sympy.polys.polynomial import Poly, SymbolsError, \
    PolynomialError, CoefficientError, MultivariatePolyError
from sympy.polys.algorithms import poly_decompose, poly_sqf, poly_div
from sympy.polys.factortools import poly_factors

from sympy.ntheory import divisors
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

    from sympy.simplify import simplify

    return [ simplify(r) for r in roots ]

def roots_cubic(f):
    """Returns a list of  roots of a cubic polynomial."""
    if isinstance(f, (tuple, list)):
        a, b, c, d = f
    else:
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
    """Returns a list of roots of a quartic polynomial.

    References:
      http://planetmath.org/encyclopedia/GaloisTheoreticDerivationOfTheQuarticFormula.html
    """
    a, b, c, d, e = f.iter_all_coeffs()

    # normalize
    t = a
    a = b/t
    b = c/t
    c = d/t
    d = e/t

    p = -2*b
    q = b**2 + a*c - 4*d
    r = c**2 + a**2*d - a*b*c

    r = roots_cubic((1, p, q, r))

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
        roots.append((alpha*zeta).expand(power_base=False))

    return roots

def roots_rational(f):
    """Returns a list of rational roots of a polynomial."""

    try:
        g = f.as_integer()[1]
    except CoefficientError:
        return []

    LC_divs = divisors(int(g.LC))
    TC_divs = divisors(int(g.TC))

    if not g(S.Zero):
        zeros = [S.Zero]
    else:
        zeros = []

    for p in LC_divs:
        for q in TC_divs:
            zero = Rational(p, q)

            if not g(zero):
                zeros.append(zero)

            if not g(-zero):
                zeros.append(-zero)

    return zeros

def roots(f, *symbols, **flags):
    """Computes symbolic roots of an univariate polynomial.

       Given an univariate  polynomial f with symbolic coefficients,
       returns a dictionary with its roots and their multiplicities.

       Only roots expressible via radicals will be returned.  To get
       a complete set of roots use RootOf class or numerical methods
       instead. By default cubic and quartic formulas are used in
       the algorithm. To disable them because of unreadable output
       set cubics=False or quartics=False respectively.

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
       {y**(1/2): 1, -y**(1/2): 1}

    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    def _update_dict(result, root, k):
        if root in result:
            result[root] += k
        else:
            result[root] = k

    def _try_decompose(f):
        """Find roots using functional decomposition. """
        factors = poly_decompose(f)
        result, g = {}, factors[0]

        for i, h in enumerate(poly_sqf(g)):
            for r in _try_heuristics(h):
                _update_dict(result, r, i+1)

        for factor in factors[1:]:
            last, result = result.copy(), {}

            for last_r, i in last.iteritems():
                g = factor.sub_term(last_r, (0,))

                for j, h in enumerate(poly_sqf(g)):
                    for r in _try_heuristics(h):
                        _update_dict(result, r, i*(j+1))

        return result

    def _try_heuristics(f):
        """Find roots using formulas and some tricks. """
        if f.length == 1:
            if f.is_constant:
                return []
            else:
                return [S(0)] * f.degree

        if f.length == 2:
            if f.degree == 1:
                return roots_linear(f)
            else:
                return roots_binomial(f)

        x, result = f.symbols[0], []

        for i in [S(-1), S(1)]:
            if f(i).expand().is_zero:
                f = poly_div(f, x-i)[0]
                result.append(i)
                break

        n = f.degree

        if n == 1:
            result += roots_linear(f)
        elif n == 2:
            result += roots_quadratic(f)
        elif n == 3 and flags.get('cubics', True):
            result += roots_cubic(f)
        elif n == 4 and flags.get('quartics', False):
            result += roots_quartic(f)

        return result

    multiple = flags.get('multiple', False)

    if f.length == 1:
        if f.is_constant:
            if multiple:
                return []
            else:
                return {}
        else:
            result = { S(0) : f.degree }
    else:
        (k,), f = f.as_reduced()

        if k == 0:
            zeros = {}
        else:
            zeros = { S(0) : k }

        result = {}

        if f.length == 2:
            if f.degree == 1:
                result[roots_linear(f)[0]] = 1
            else:
                for r in roots_binomial(f):
                    _update_dict(result, r, 1)
        elif f.degree == 2:
            for r in roots_quadratic(f):
                _update_dict(result, r, 1)
        else:
            try:
                _, factors = poly_factors(f)

                if len(factors) == 1 and factors[0][1] == 1:
                    raise CoefficientError

                for factor, k in factors:
                    for r in _try_heuristics(factor):
                        _update_dict(result, r, k)
            except CoefficientError:
                result = _try_decompose(f)

        result.update(zeros)

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
            raise ValueError("Invalid domain: %s" % domain)

        for zero in dict(result).iterkeys():
            if not query(zero):
                del result[zero]

    predicate = flags.get('predicate', None)

    if predicate is not None:
        for zero in dict(result).iterkeys():
            if not predicate(zero):
                del result[zero]

    if not multiple:
        return result
    else:
        zeros = []

        for zero, k in result.iteritems():
            zeros.extend([zero] * k)

        return zeros

def poly_root_factors(f, *symbols, **flags):
    """Returns all factors of an univariate polynomial.

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> factors = poly_root_factors(x**2-y, x)

       >>> set(f.as_basic() for f in factors)
       set([x + y**(1/2), x - y**(1/2)])
    """
    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    if f.is_multivariate:
        raise MultivariatePolyError(f)
    else:
        x = f.symbols[0]

    if 'multiple' in flags:
        del flags['multiple']

    zeros = roots(f, **flags)

    if not zeros:
        return [f]
    else:
        factors, N = [], 0

        for r, n in zeros.iteritems():
            h = Poly([(S.One, 1), (-r, 0)], x)
            factors, N = factors + [h]*n, N + n

        if N < f.degree:
            g = reduce(lambda p,q: p*q, factors)
            factors.append(poly_div(f, g)[0])

        return factors

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
        raise SymbolsError("Redundant symbols were given")

    if f.is_multivariate:
        raise MultivariatePolyError(f)
    else:
        f = f.as_squarefree()

    sturm = [f, f.diff()]

    while not sturm[-1].is_zero:
        sturm.append(-poly_div(sturm[-2], sturm[-1])[1])

    return sturm[:-1]

def number_of_real_roots(f, *symbols, **flags):
    """Returns the number of distinct real roots of f in (inf, sup].

       >>> from sympy import *
       >>> x,y = symbols('xy')

       >>> f = Poly(x**2 - 1, x)

       Count real roots in the (-oo, oo) interval:

       >>> number_of_real_roots(f)
       2

       Count real roots in the (0, 2) interval:

       >>> number_of_real_roots(f, inf=0, sup=2)
       1

       Count real roots in the (sqrt(2), oo) interval:

       >>> number_of_real_roots(f, inf=sqrt(2))
       0

       For more information on the implemented algorithm refer to:

       [1] J.H. Davenport, Y. Siret, E. Tournier, Computer Algebra
           Systems and Algorithms for Algebraic Computation,
           Academic Press, London, 1988, pp. 124-128
    """

    def sign_changes(seq):
        count = 0

        for i in xrange(1, len(seq)):
            if (seq[i-1] < 0 and seq[i] >= 0) or \
               (seq[i-1] > 0 and seq[i] <= 0):
                count += 1

        return count

    if not isinstance(f, Poly):
        f = Poly(f, *symbols)
    elif symbols:
        raise SymbolsError("Redundant symbols were given")

    if f.is_multivariate:
        raise MultivariatePolyError(f)

    if f.degree < 1:
        return 0

    inf = flags.get('inf', None)

    if inf is not None:
        inf = sympify(inf)

        if not inf.is_number:
            raise ValueError("Not a number: %s" % inf)
        elif abs(inf) is S.Infinity:
            inf = None

    sup = flags.get('sup', None)

    if sup is not None:
        sup = sympify(sup)

        if not sup.is_number:
            raise ValueError("Not a number: %s" % sup)
        elif abs(sup) is S.Infinity:
            sup = None

    sturm = poly_sturm(f)

    if inf is None:
        signs_inf = sign_changes([ s.LC * (-1)**s.LM[0] for s in sturm ])
    else:
        signs_inf = sign_changes([ s(inf) for s in sturm ])

    if sup is None:
        signs_sup = sign_changes([ s.LC for s in sturm ])
    else:
        signs_sup = sign_changes([ s(sup) for s in sturm ])

    return abs(signs_inf - signs_sup)

_exact_roots_cache = {}

def _exact_roots(f):
    if f in _exact_roots_cache:
        zeros = _exact_roots_cache[f]
    else:
        exact, zeros = roots(f), []

        for zero, k in exact.iteritems():
            zeros += [zero] * k

        _exact_roots_cache[f] = zeros

    return zeros

class RootOf(Basic):
    """Represents n-th root of an univariate polynomial. """

    def __new__(cls, f, index):
        if isinstance(f, RootsOf):
            f = f.poly
        elif not isinstance(f, Poly):
            raise PolynomialError("%s is not a polynomial" % f)

        if f.is_multivariate:
            raise MultivariatePolyError(f)

        if index < 0 or index >= f.degree:
            raise IndexError("Index must be in [0, %d] range" % (f.degree-1))
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

class RootsOf(Basic):
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
            raise SymbolsError("Redundant symbols were given")

        if f.is_multivariate:
            raise MultivariatePolyError(f)

        return Basic.__new__(cls, f)

    @property
    def poly(self):
        return self._args[0]

    @property
    def count(self):
        return self.poly.degree

    def roots(self):
        """Iterates over all roots: exact and formal. """
        exact = _exact_roots(self.poly)

        for root in exact:
            yield root

        for j in range(len(exact), self.count):
            yield RootOf(self.poly, j)

    def exact_roots(self):
        """Iterates over exact roots only. """
        exact = _exact_roots(self.poly)

        for root in exact:
            yield root

    def formal_roots(self):
        """Iterates over formal roots only. """
        exact = _exact_roots(self.poly)

        for j in range(len(exact), self.count):
            yield RootOf(self.poly, j)

    def __call__(self, index):
        return RootOf(self.poly, index)

    def atoms(self, *args, **kwargs):
        return self.poly.atoms(*args, **kwargs)


class RootSum(Basic):
    """Represents a sum of all roots of an univariate polynomial. """

    def __new__(cls, f, *args, **flags):
        if not hasattr(f, '__call__'):
            raise TypeError("%s is not a callable object" % f)

        roots = RootsOf(*args)

        if not flags.get('evaluate', True):
            return Basic.__new__(cls, f, roots)
        else:
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

            for root in self.roots.roots():
                result += self.function(root)

            return result

