
from sympy import Basic, S, C, Symbol, Mul, I, pi, exp, sqrt
from sympy.core import sympify
from sympy.core.numbers import Rational, Integer, gcd
from sympy.core.methods import NoRelMeths, ArithMeths
from sympy.polynomials import quo, factor_, PolynomialException
from sympy.numerics.optimize import polyroots

def linear(a, b):
    return [-b/a]

def quadratic(a, b, c):
    d = (b**2 - 4*a*c)**S.Half

    R = [(-b + d) / (2*a),
         (-b - d) / (2*a)]

    return [ r.expand() for r in R ]

def cubic(a, b, c, d):
    p = (3*c/a - (b/a)**2) / 9
    q = (9*b*c/a**2 - 27*d/a - 2*(b/a)**3) / 54

    d = (p**3 + q**2)**S.Half

    s = (q - d)**Rational(1,3)
    t = (q - d)**Rational(1,3)

    I, J = s + t, s - t

    u, v = b/(3*a), I*sqrt(3)/2

    R = [I - u,
         -I/2 - u + v*J,
         -I/2 - u - v*J]

    return [ r.expand() for r in R ]

def quartic(a, b, c, d, e):
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

    R = [(C + E) / 2, (C - E) / 2,
         (D + F) / 2, (D - F) / 2]

    return [ r.expand() for r in R ]

def arbitrary(poly):
    g = reduce(gcd, [ int(e) for _, e in poly.coeffs ])

    if g in [0, 1]:
        return []

    n = int(poly.degree() / g)

    if n not in [1, 2, 3]:
        return []

    poly = C.Polynomial(coeffs=[ (c, e/g) for c, e in poly.coeffs ],
        var=poly.var, order=poly.order)

    roots_of_unity = [ exp(2*i*pi*I/g) for i in range(g) ]

    solutions = []

    for zeta in roots_of_unity:
        zeta = zeta.expand(complex=True)

        for r in roots(poly):
            R = r**Rational(1,g)

            if not r.atoms(Symbol):
                R = R.expand(complex=True)

            solutions.append((zeta * R).expand())

    return solutions

def roots(poly, x=None, domain=None, **kwargs):
    if not isinstance(poly, C.Polynomial):
        try:
            poly = poly.as_polynomial(x)
        except PolynomialException:
            return []

    n = int(poly.degree())

    if n == 0:
        return []

    if domain == 'F':
        if len(poly.atoms(Symbol)) == 1:
            return sorted([ sympify(complex(r)) \
                for r in polyroots(poly.as_basic())[0] ])
        else:
            raise ValueError("'F' domain requires numeric coefficients")

    roots = arbitrary(poly)

    if roots == []:
        factors = [ poly ]

        if n < 5 or kwargs.get('factor', False):
            try:
                factors = factor_.factor(poly)
            except PolynomialException:
                pass

        handlers = {
            1 : linear,
            2 : quadratic,
            3 : cubic,
            4 : quartic
        }

        for factor in factors:
            n = int(factor.degree())

            if n == 0:
                continue

            if n in handlers:
                coeffs = [ factor.nth_coeff(i) for i in range(n+1) ]
                roots += handlers[n](*reversed(coeffs))
            else:
                roots += arbitrary(factor)

    Zr, Qr, Rr, Ir, Cr = [], [], [], [], []

    for root in roots:
        if root.is_Integer:
            Zr.append(root)
        elif root.is_Rational:
            Qr.append(root)
        elif root.is_real:
            Rr.append(root)
        elif root.is_imaginary:
            Ir.append(root)
        else:
            Cr.append(root)

    Zr, Qr, Rr = map(sorted, (Zr, Qr, Rr))

    domains = {
        'Z' : Zr,
        'Q' : Zr+Qr,
        'R' : Zr+Qr+Rr,
        'I' : Ir,
        'C' : Zr+Qr+Rr+Ir+Cr
    }

    if domain is None:
        roots = domains['C']
    else:
        try:
            roots = domains[domain]
        except KeyError:
            raise ValueError("Invalid domain: %s" % domain)

    return roots

def factors(poly, x=None, domain=None, **kwargs):
    _roots = roots(poly, x, domain, **kwargs)

    if isinstance(poly, C.Polynomial):
        poly = poly.as_basic()

    if not _roots:
        return [ poly ]
    else:
        factored = [ x-r for r in _roots ]
        f = quo(poly, Mul(*factored), x)

        if f.has(x):
            return factored + [ f ]
        else:
            return factored

_exact_roots_cache = {}
_float_roots_cache = {}

def _exact_roots(poly):
    try:
        exact = _exact_roots_cache[poly]
    except KeyError:
        exact = roots(poly, domain='C')
        _exact_roots_cache[poly] = exact

    return exact

def _float_roots(poly):
    try:
        approx = _float_roots_cache[poly]
    except KeyError:
        approx = roots(poly, domain='F')
        _float_roots_cache[poly] = approx

    return approx

class RootOf(Basic, NoRelMeths, ArithMeths):
    """Represents all roots or a specific root of a polynomial.

       RootOf class is dual to support handling of all roots
       of a given polynomial or just a single root specified
       using integer index. In both cases roots can be computed
       exactly, if possible, or using numerical approximation.

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> prec = Basic.set_precision(6)

       >>> RootOf(2 + x**2, x).count()
       2

       >>> RootOf(2 + x**2, x, 0)
       I*2**(1/2)
       >>> RootOf(2 + x**2, x, 1)
       -I*2**(1/2)

       >>> list(RootOf(2 + x**2, x).roots())
       [I*2**(1/2), -I*2**(1/2)]

       >>> roots = list(RootOf(1 + x + x**7, x).roots())
       >>> print roots[0]
       RootOf(1 + x + x**7, x, 0)
       >>> print roots[1]
       RootOf(1 + x + x**7, x, 1)

       and so on ...

       >>> RootOf(1 + x + x**7, x, 1).evalf()
       -0.705298 - 0.637624*I

       >>> _ = Basic.set_precision(prec)

    """

    def __new__(cls, poly, x, index=None, **assumptions):
        try:
            poly = poly.as_polynomial(x)
        except PolynomialException:
            raise TypeError("'RootOf' operates on polynomials")

        if index is not None:
            if index < 0 or index >= poly.degree():
                raise IndexError("Index must be in "
                    "[0, %d] range" % int(poly.degree()-1))
            else:
                exact = _exact_roots(poly)

                if index < len(exact):
                    return exact[index]

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (poly, x, index)

        return obj

    @property
    def poly(self):
        return self._args[0]

    @property
    def x(self):
        return self._args[1]

    @property
    def index(self):
        return self._args[2]

    @property
    def represents_all(self):
        return self.index is None

    @property
    def represents_one(self):
        return self.index is not None

    def tostr(self, level=0):
        poly, x = (self.poly.tostr(), self.x.tostr())

        if self.represents_one:
            return "RootOf(%s, %s, %d)" % (poly, x, self.index)
        else:
            return "RootOf(%s, %s)" % (poly, x)

    def count(self):
        if self.represents_one:
            return 1
        else:
            return int(self.poly.degree())

    def roots(self):
        if self.represents_one:
            yield self
        else:
            exact = _exact_roots(self.poly)

            for i, root in enumerate(exact):
                yield root

            for i in range(len(exact), self.count()):
                yield RootOf(self.poly, self.x, i)

    def _eval_evalf(self):
        if self.represents_one:
            approx = _float_roots(self.poly)
            return approx[self.index]
        else:
            raise TypeError("Use 'approx' method to obtain all roots")

    def approx(self):
        if self.represents_one:
            yield self.evalf()
        else:
            approx = _float_roots(self.poly)

            for root in approx:
                yield root

    def atoms(self, *args, **kwargs):
        return self.poly.atoms(*args, **kwargs)
