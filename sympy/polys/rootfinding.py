from sympy.core.symbol import Symbol
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
from sympy.functions import exp, sqrt

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

    return roots

def roots_cubic(f):
    """Returns a list of roots of a cubic polynomial."""
    one, a, b, c = f.as_monic().iter_all_coeffs()

    if c is S.Zero:
        x1, x2 = roots([1,a,b], multiple = True)
        return [x1, S.Zero, x2]

    p = b - a**2/3
    q = c - a*b/3 + 2*a**3/27

    pon3 = p/3
    aon3 = a/3

    if p is S.Zero:
        if q is S.Zero:
            return [-aon3] * 3
        else:
            u1 = q**Rational(1, 3)
    elif q is S.Zero:
        y1, y2 = roots([1, 0, p], multiple=True)
        return [tmp - aon3 for tmp in [y1, S.Zero, y2]]
    else:
        u1 = (q/2 + (q**2/4 + pon3**3)**S.Half)**Rational(1, 3)

    coeff = S.ImaginaryUnit*3**S.Half / 2

    u2 = u1*(-S.Half + coeff)
    u3 = u1*(-S.Half - coeff)

    soln = [
        -u1 + pon3/u1 - aon3,
        -u2 + pon3/u2 - aon3,
        -u3 + pon3/u3 - aon3
    ]

    return soln

def roots_quartic(f):
    """Returns a list of roots of a quartic polynomial.

    There are many references for solving quartic expressions available [1-5].
    This reviewer has found that many of them require one to select from among
    2 or more possible sets of solutions and that some solutions work when one
    is searching for real roots but don't work when searching for complex roots
    (though this is not always stated clearly). The following routine has been
    tested and found to be correct for 0, 2 or 4 complex roots.

    The quasisymmetric case solution[6] looks for quartics that have the form
    x**4 + A*x**3 + B*x**2 + C*x + D = 0 where (C/A)**2 = D.

    NOTE: There is not a single symbolic solution that is valid for all
    possible values of A, B, C and D. There are 4 sets of solutions possible
    based on the code below. These solutions (determined by the values of the
    reduced quartic coefficients, a = B/A, b = C/A, c = D/A, and d = E/A, are:
      1) f = c + a * (a**2 / 8 - b / 2) == 0
      2) g = d - a * (a * (3 * a**2 / 256 - b / 16) + c / 4) = 0
      3) if f != 0 and g !=0 but
        a) p = -d + a*c/4 - b**2/12 = 0 and
           q = b*d/3 + a*b*c/24 - c**2/8 - d*a**2/8 - b**3/108 >= 0
           then u (see code) will be zero, otherwise
        b) u != 0

            Schematically, it looks like this:

                f = 0    g = 0           f != 0 and g != 0
                  |        |                   /   \
                  |        |       p=0 and q>=0     p!=0 or q<0
                  |        |       (i.e. u = 0     (i.e. u != 0)
                  |        |            |                |
                4 solns  4 solns      4 solns          4 solns
                                             (default symbolic solution)

    These branches do not count the special cases that have a particularly
    simple form of the roots. Those special cases can often be solved by the
    general procedure. In the test suite, there are tests for those cases,
    and each of them marked with "general soln ok, too" in the code below
    still gave the same for the tests.

    Example::
    >>> from sympy import var, Poly
    >>> from sympy.polys.rootfinding import roots_quartic
    >>> x = var('x')
    >>> r = roots_quartic(Poly(x**4 -6*x**3 +17*x**2 -26*x +20, x))
    >>> # 4 complex roots: 1+-I*sqrt(3), 2+-I
    >>> sorted(str(tmp.evalf(n=2)) for tmp in r)
    ['1.0 + 1.7*I', '1.0 - 1.7*I', '2.0 + I', '2.0 - 1.0*I']

    References:
      [1] http://mathforum.org/dr.math/faq/faq.cubic.equations.html
      [2] http://en.wikipedia.org/wiki/Quartic_function#Summary_of_Ferrari.27s_method
      [3] http://planetmath.org/encyclopedia/GaloisTheoreticDerivationOfTheQuarticFormula.html
      [4] http://staff.bath.ac.uk/masjhd/JHD-CA.pdf
      [5] http://www.albmath.org/files/Math_5713.pdf
      [6] http://www.statemaster.com/encyclopedia/Quartic-equation
          #Other_particular_case:_Quasi-symmetric_equations

    """

    # normalized coefficients
    one, a, b, c, d = f.as_monic().iter_all_coeffs()

    if d is S.Zero: # general soln ok, too
        return [S.Zero] + roots([1, a, b, c], multiple = True)
    elif (c/a)**2 == d: # general solution ok, too
        m = sqrt(d)
        z = Symbol('z', dummy=True)
        x = f.symbols[0]
        z1, z2 = roots_quadratic(Poly(z**2+a*z+b-2*m, z))
        f = Poly(x**2 - z*x + m, x)
        return roots_quadratic(f.subs(z, z1)) + \
               roots_quadratic(f.subs(z, z2))
    else:
        a2 = a ** 2
        e = b - 3 * a2 / 8
        f = c + a * (a2 / 8 - b / 2)
        g = d - a * (a * (3 * a2 / 256 - b / 16) + c / 4)
        aon4 = a / 4
        ans = []

        if f.is_zero: # general solution not valid if f = 0.
            y1, y2 = [tmp ** S.Half for tmp in
                      roots([1, e, g], multiple = True)]
            return [tmp - aon4 for tmp in [-y1, -y2, y1, y2]]
        if g.is_zero: # general solution not valid if g = 0
            y = [S.Zero] + roots([1, 0, e, f], multiple = True)
            return [tmp - aon4 for tmp in y]
        else:
            TH = Rational(1, 3)
            p = -e**2/12 - g
            q = -e**3/108 + e*g/3 - f**2/8
            root = sqrt(q**2/4 + p**3/27)
            rr = [-q/2 + s*root for s in [1]] # in [1,-1], either will do, so pick 1
            for r in rr:
                uu = [r**TH] # solve(x**3-r,x), any one will do, so take primary
                for u in uu:
                    if u.is_zero:
                        y = -5*e/6 + u - q**TH
                    else:
                        y = -5*e/6 + u - p/u/3
                    w = sqrt(e + 2*y)
                    # try to return the real solutions first if the
                    # sign can be determined. The term tested is the
                    # argument to `root` below.
                    arg1 = 3*e + 2*y
                    arg2 = 2*f/w
                    if -(arg1 + arg2) > S.Zero:
                        ss = [-1, 1]
                    else:
                        ss= [1, -1]
                    for s in ss:
                        root = sqrt(-(arg1 + s*arg2))
                        # return the more negative root first;
                        # note that t and s have opposite signs
                        # in the formula: s=+/-1 while t=-/+1.
                        if aon4 > S.Zero:
                            tt = [-1, 1]
                        else:
                            tt = [1, -1]
                        for t in tt:
                            ans.append((s*w - t*root)/2 - aon4)
    return ans

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
    """Computes symbolic roots of a univariate polynomial.

       Given a univariate polynomial f with symbolic coefficients (or
       a list of the polynomial's coefficients), returns a dictionary
       with its roots and their multiplicities.

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

       >>> from sympy import Poly, roots
       >>> from sympy.abc import x, y

       >>> roots(x**2 - 1, x)
       {1: 1, -1: 1}

       >>> p = Poly(x**2-1, x)
       >>> roots(p)
       {1: 1, -1: 1}

       >>> m = Poly(x**2-y, x, y) # multivariate
       >>> roots(Poly(m, x)) # same expr, but now univariate
       {y**(1/2): 1, -y**(1/2): 1}

       >>> roots([1, 0, -1]) # LIST of coefficients
       {1: 1, -1: 1}

       >>> roots(x**2 - y, x)
       {y**(1/2): 1, -y**(1/2): 1}

    """

    if isinstance(f, list):
        if symbols:
            raise SymbolsError("Redundant symbol(s) given.")
        f = Poly(f, Symbol('x', dummy=True))

    if not isinstance(f, Poly):
        if len(symbols) > 1:
            raise SymbolsError("Too many symbols given.")
        f = Poly(f, *(symbols or S(f).atoms(Symbol)))

    if len(symbols) > 1:
        raise SymbolsError("Too many symbols given.")

    if f.is_multivariate:
        if len(symbols) == 1:
            f = Poly(f, *symbols)
        else:
            raise MultivariatePolyError(f)

    assert f.is_univariate

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
        elif n == 4 and flags.get('quartics', True):
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
    """Returns all factors of a univariate polynomial.

       >>> from sympy.polys.rootfinding import poly_root_factors
       >>> from sympy.abc import x, y

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

       Given a univariate, square-free polynomial f(x) returns an
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

       >>> from sympy.polys.rootfinding import number_of_real_roots
       >>> from sympy.abc import x, y
       >>> from sympy import Poly, sqrt

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
    """Represents n-th root of a univariate polynomial. """

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
    """Represents all roots of a univariate polynomial.

       >>> from sympy import roots, RootsOf
       >>> from sympy.abc import x, y

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
    """Represents a sum of all roots of a univariate polynomial. """

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
