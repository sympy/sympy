"""Implementation of RootOf class and related tools. """

from __future__ import print_function, division

from sympy.core import (S, Expr, Integer, Float, I, Add, Lambda, symbols,
        sympify, Rational)

from sympy.polys.polytools import Poly, PurePoly, factor
from sympy.polys.rationaltools import together
from sympy.polys.polyfuncs import symmetrize, viete

from sympy.polys.rootisolation import (
    dup_isolate_complex_roots_sqf,
    dup_isolate_real_roots_sqf)

from sympy.polys.polyroots import (
    roots_linear, roots_quadratic, roots_binomial,
    preprocess_roots, roots)

from sympy.polys.polyerrors import (
    MultivariatePolynomialError,
    GeneratorsNeeded,
    PolynomialError,
    DomainError)

from sympy.polys.domains import QQ

from sympy.mpmath import mp, mpf, mpc, findroot
from sympy.mpmath.libmp.libmpf import prec_to_dps

from sympy.utilities import lambdify, public

from sympy.core.compatibility import xrange

_reals_cache = {}
_complexes_cache = {}

@public
class RootOf(Expr):
    """Represents ``k``-th root of a univariate polynomial. """

    __slots__ = ['poly', 'index']
    is_complex = True

    def __new__(cls, f, x, index=None, radicals=True, expand=True):
        """Construct a new ``RootOf`` object for ``k``-th root of ``f``. """
        x = sympify(x)

        if index is None and x.is_Integer:
            x, index = None, x
        else:
            index = sympify(index)

        if index.is_Integer:
            index = int(index)
        else:
            raise ValueError("expected an integer root index, got %d" % index)

        poly = PurePoly(f, x, greedy=False, expand=expand)

        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are allowed")

        degree = poly.degree()

        if degree <= 0:
            raise PolynomialError("can't construct RootOf object for %s" % f)

        if index < -degree or index >= degree:
            raise IndexError("root index out of [%d, %d] range, got %d" %
                             (-degree, degree - 1, index))
        elif index < 0:
            index += degree

        dom = poly.get_domain()

        if not dom.is_Exact:
            poly = poly.to_exact()

        roots = cls._roots_trivial(poly, radicals)

        if roots is not None:
            return roots[index]

        coeff, poly = preprocess_roots(poly)
        dom = poly.get_domain()

        if not dom.is_ZZ:
            raise NotImplementedError("RootOf is not supported over %s" % dom)

        root = cls._indexed_root(poly, index)
        return coeff*cls._postprocess_root(root, radicals)

    @classmethod
    def _new(cls, poly, index):
        """Construct new ``RootOf`` object from raw data. """
        obj = Expr.__new__(cls)

        obj.poly = poly
        obj.index = index

        return obj

    def _hashable_content(self):
        return (self.poly, self.index)

    @property
    def expr(self):
        return self.poly.as_expr()

    @property
    def args(self):
        return (self.expr, Integer(self.index))

    @property
    def free_symbols(self):
        # RootOf currently only works with univariate expressions and although
        # the poly attribute is often a PurePoly, sometimes it is a Poly. In
        # either case no free symbols should be reported.
        return set()

    def _eval_is_real(self):
        """Return ``True`` if the root is real. """
        return self.index < len(_reals_cache[self.poly])

    @classmethod
    def real_roots(cls, poly, radicals=True):
        """Get real roots of a polynomial. """
        return cls._get_roots("_real_roots", poly, radicals)

    @classmethod
    def all_roots(cls, poly, radicals=True):
        """Get real and complex roots of a polynomial. """
        return cls._get_roots("_all_roots", poly, radicals)

    @classmethod
    def _get_reals_sqf(cls, factor):
        """Compute real root isolating intervals for a square-free polynomial. """
        if factor in _reals_cache:
            real_part = _reals_cache[factor]
        else:
            _reals_cache[factor] = real_part = \
                dup_isolate_real_roots_sqf(
                    factor.rep.rep, factor.rep.dom, blackbox=True)

        return real_part

    @classmethod
    def _get_complexes_sqf(cls, factor):
        """Compute complex root isolating intervals for a square-free polynomial. """
        if factor in _complexes_cache:
            complex_part = _complexes_cache[factor]
        else:
            _complexes_cache[factor] = complex_part = \
                dup_isolate_complex_roots_sqf(
                    factor.rep.rep, factor.rep.dom, blackbox=True)

        return complex_part

    @classmethod
    def _get_reals(cls, factors):
        """Compute real root isolating intervals for a list of factors. """
        reals = []

        for factor, k in factors:
            real_part = cls._get_reals_sqf(factor)
            reals.extend([ (root, factor, k) for root in real_part ])

        return reals

    @classmethod
    def _get_complexes(cls, factors):
        """Compute complex root isolating intervals for a list of factors. """
        complexes = []

        for factor, k in factors:
            complex_part = cls._get_complexes_sqf(factor)
            complexes.extend([ (root, factor, k) for root in complex_part ])

        return complexes

    @classmethod
    def _reals_sorted(cls, reals):
        """Make real isolating intervals disjoint and sort roots. """
        cache = {}

        for i, (u, f, k) in enumerate(reals):
            for j, (v, g, m) in enumerate(reals[i + 1:]):
                u, v = u.refine_disjoint(v)
                reals[i + j + 1] = (v, g, m)

            reals[i] = (u, f, k)

        reals = sorted(reals, key=lambda r: r[0].a)

        for root, factor, _ in reals:
            if factor in cache:
                cache[factor].append(root)
            else:
                cache[factor] = [root]

        for factor, roots in cache.items():
            _reals_cache[factor] = roots

        return reals

    @classmethod
    def _complexes_sorted(cls, complexes):
        """Make complex isolating intervals disjoint and sort roots. """
        cache = {}

        for i, (u, f, k) in enumerate(complexes):
            for j, (v, g, m) in enumerate(complexes[i + 1:]):
                u, v = u.refine_disjoint(v)
                complexes[i + j + 1] = (v, g, m)

            complexes[i] = (u, f, k)

        complexes = sorted(complexes, key=lambda r: (r[0].ax, r[0].ay))

        for root, factor, _ in complexes:
            if factor in cache:
                cache[factor].append(root)
            else:
                cache[factor] = [root]

        for factor, roots in cache.items():
            _complexes_cache[factor] = roots

        return complexes

    @classmethod
    def _reals_index(cls, reals, index):
        """Map initial real root index to an index in a factor where the root belongs. """
        i = 0

        for j, (_, factor, k) in enumerate(reals):
            if index < i + k:
                poly, index = factor, 0

                for _, factor, _ in reals[:j]:
                    if factor == poly:
                        index += 1

                return poly, index
            else:
                i += k

    @classmethod
    def _complexes_index(cls, complexes, index):
        """Map initial complex root index to an index in a factor where the root belongs. """
        index, i = index, 0

        for j, (_, factor, k) in enumerate(complexes):
            if index < i + k:
                poly, index = factor, 0

                for _, factor, _ in complexes[:j]:
                    if factor == poly:
                        index += 1

                index += len(_reals_cache[poly])

                return poly, index
            else:
                i += k

    @classmethod
    def _count_roots(cls, roots):
        """Count the number of real or complex roots including multiplicites. """
        return sum([ k for _, _, k in roots ])

    @classmethod
    def _indexed_root(cls, poly, index):
        """Get a root of a composite polynomial by index. """
        (_, factors) = poly.factor_list()

        reals = cls._get_reals(factors)
        reals_count = cls._count_roots(reals)

        if index < reals_count:
            reals = cls._reals_sorted(reals)
            return cls._reals_index(reals, index)
        else:
            complexes = cls._get_complexes(factors)
            complexes = cls._complexes_sorted(complexes)
            return cls._complexes_index(complexes, index - reals_count)

    @classmethod
    def _real_roots(cls, poly):
        """Get real roots of a composite polynomial. """
        (_, factors) = poly.factor_list()

        reals = cls._get_reals(factors)
        reals = cls._reals_sorted(reals)
        reals_count = cls._count_roots(reals)

        roots = []

        for index in xrange(0, reals_count):
            roots.append(cls._reals_index(reals, index))

        return roots

    @classmethod
    def _all_roots(cls, poly):
        """Get real and complex roots of a composite polynomial. """
        (_, factors) = poly.factor_list()

        reals = cls._get_reals(factors)
        reals = cls._reals_sorted(reals)
        reals_count = cls._count_roots(reals)

        roots = []

        for index in xrange(0, reals_count):
            roots.append(cls._reals_index(reals, index))

        complexes = cls._get_complexes(factors)
        complexes = cls._complexes_sorted(complexes)
        complexes_count = cls._count_roots(complexes)

        for index in xrange(0, complexes_count):
            roots.append(cls._complexes_index(complexes, index))

        return roots

    @classmethod
    def _roots_trivial(cls, poly, radicals):
        """Compute roots in linear, quadratic and binomial cases. """
        if poly.degree() == 1:
            return roots_linear(poly)

        if not radicals:
            return None

        if poly.degree() == 2:
            return roots_quadratic(poly)
        elif poly.length() == 2 and poly.TC():
            return roots_binomial(poly)
        else:
            return None

    @classmethod
    def _preprocess_roots(cls, poly):
        """Take heroic measures to make ``poly`` compatible with ``RootOf``. """
        dom = poly.get_domain()

        if not dom.is_Exact:
            poly = poly.to_exact()

        coeff, poly = preprocess_roots(poly)
        dom = poly.get_domain()

        if not dom.is_ZZ:
            raise NotImplementedError("RootOf is not supported over %s" % dom)

        return coeff, poly

    @classmethod
    def _postprocess_root(cls, root, radicals):
        """Return the root if it is trivial or a ``RootOf`` object. """
        poly, index = root
        roots = cls._roots_trivial(poly, radicals)

        if roots is not None:
            return roots[index]
        else:
            return cls._new(poly, index)

    @classmethod
    def _get_roots(cls, method, poly, radicals):
        """Return postprocessed roots of specified kind. """
        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are allowed")

        coeff, poly = cls._preprocess_roots(poly)
        roots = []

        for root in getattr(cls, method)(poly):
            roots.append(coeff*cls._postprocess_root(root, radicals))

        return roots

    def _get_interval(self):
        """Internal function for retrieving isolation interval from cache. """
        if self.is_real:
            return _reals_cache[self.poly][self.index]
        else:
            reals_count = len(_reals_cache[self.poly])
            return _complexes_cache[self.poly][self.index - reals_count]

    def _set_interval(self, interval):
        """Internal function for updating isolation interval in cache. """
        if self.is_real:
            _reals_cache[self.poly][self.index] = interval
        else:
            reals_count = len(_reals_cache[self.poly])
            _complexes_cache[self.poly][self.index - reals_count] = interval

    def _eval_evalf(self, prec):
        """Evaluate this complex root to the given precision. """
        _prec, mp.prec = mp.prec, prec

        try:
            func = lambdify(self.poly.gen, self.expr)

            interval = self._get_interval()
            if not self.is_real:
                # For complex intervals, we need to keep refining until the
                # imaginary interval is disjunct with other roots, that is,
                # until both ends get refined.
                ay = interval.ay
                by = interval.by
                while interval.ay == ay or interval.by == by:
                    interval = interval.refine()

            while True:
                if self.is_real:
                    x0 = mpf(str(interval.center))
                else:
                    x0 = mpc(*map(str, interval.center))

                try:
                    root = findroot(func, x0)
                    # If the (real or complex) root is not in the 'interval',
                    # then keep refining the interval. This happens if findroot
                    # accidentally finds a different root outside of this
                    # interval because our initial estimate 'x0' was not close
                    # enough.
                    if self.is_real:
                        a = mpf(str(interval.a))
                        b = mpf(str(interval.b))
                        # This is needed due to the bug #3364:
                        a, b = min(a, b), max(a, b)
                        if not (a < root < b):
                            raise ValueError("Root not in the interval.")
                    else:
                        ax = mpf(str(interval.ax))
                        bx = mpf(str(interval.bx))
                        ay = mpf(str(interval.ay))
                        by = mpf(str(interval.by))
                        # This is needed due to the bug #3364:
                        ax, bx = min(ax, bx), max(ax, bx)
                        ay, by = min(ay, by), max(ay, by)
                        if not (ax < root.real < bx and ay < root.imag < by):
                            raise ValueError("Root not in the interval.")
                except ValueError:
                    interval = interval.refine()
                    continue
                else:
                    break
        finally:
            mp.prec = _prec

        return Float._new(root.real._mpf_, prec) + I*Float._new(root.imag._mpf_, prec)

    def eval_rational(self, tol):
        """
        Returns a Rational approximation to ``self`` with the tolerance ``tol``.

        This method uses bisection, which is very robust and it will always
        converge. The returned Rational instance will be at most 'tol' from the
        exact root.

        The following example first obtains Rational approximation to 1e-7
        accuracy for all roots of the 4-th order Legendre polynomial, and then
        evaluates it to 5 decimal digits (so all digits will be correct
        including rounding):

        >>> from sympy import S, legendre_poly, Symbol
        >>> x = Symbol("x")
        >>> p = legendre_poly(4, x, polys=True)
        >>> roots = [r.eval_rational(S(1)/10**7) for r in p.real_roots()]
        >>> roots = [str(r.n(5)) for r in roots]
        >>> roots
        ['-0.86114', '-0.33998', '0.33998', '0.86114']

        """

        if not self.is_real:
            raise NotImplementedError("eval_rational() only works for real polynomials so far")
        func = lambdify(self.poly.gen, self.expr)
        interval = self._get_interval()
        a = Rational(str(interval.a))
        b = Rational(str(interval.b))
        # This is needed due to the bug #3364:
        a, b = min(a, b), max(a, b)
        return bisect(func, a, b, tol)

@public
class RootSum(Expr):
    """Represents a sum of all roots of a univariate polynomial. """

    __slots__ = ['poly', 'fun', 'auto']

    def __new__(cls, expr, func=None, x=None, auto=True, quadratic=False):
        """Construct a new ``RootSum`` instance carrying all roots of a polynomial. """
        coeff, poly = cls._transform(expr, x)

        if not poly.is_univariate:
            raise MultivariatePolynomialError(
                "only univariate polynomials are allowed")

        if func is None:
            func = Lambda(poly.gen, poly.gen)
        else:
            try:
                is_func = func.is_Function
            except AttributeError:
                is_func = False

            if is_func and 1 in func.nargs:
                if not isinstance(func, Lambda):
                    func = Lambda(poly.gen, func(poly.gen))
            else:
                raise ValueError(
                    "expected a univariate function, got %s" % func)

        var, expr = func.variables[0], func.expr

        if coeff is not S.One:
            expr = expr.subs(var, coeff*var)

        deg = poly.degree()

        if not expr.has(var):
            return deg*expr

        if expr.is_Add:
            add_const, expr = expr.as_independent(var)
        else:
            add_const = S.Zero

        if expr.is_Mul:
            mul_const, expr = expr.as_independent(var)
        else:
            mul_const = S.One

        func = Lambda(var, expr)

        rational = cls._is_func_rational(poly, func)
        (_, factors), terms = poly.factor_list(), []

        for poly, k in factors:
            if poly.is_linear:
                term = func(roots_linear(poly)[0])
            elif quadratic and poly.is_quadratic:
                term = sum(map(func, roots_quadratic(poly)))
            else:
                if not rational or not auto:
                    term = cls._new(poly, func, auto)
                else:
                    term = cls._rational_case(poly, func)

            terms.append(k*term)

        return mul_const*Add(*terms) + deg*add_const

    @classmethod
    def _new(cls, poly, func, auto=True):
        """Construct new raw ``RootSum`` instance. """
        obj = Expr.__new__(cls)

        obj.poly = poly
        obj.fun = func
        obj.auto = auto

        return obj

    @classmethod
    def new(cls, poly, func, auto=True):
        """Construct new ``RootSum`` instance. """
        if not func.expr.has(*func.variables):
            return func.expr

        rational = cls._is_func_rational(poly, func)

        if not rational or not auto:
            return cls._new(poly, func, auto)
        else:
            return cls._rational_case(poly, func)

    @classmethod
    def _transform(cls, expr, x):
        """Transform an expression to a polynomial. """
        poly = PurePoly(expr, x, greedy=False)
        return preprocess_roots(poly)

    @classmethod
    def _is_func_rational(cls, poly, func):
        """Check if a lambda is areational function. """
        var, expr = func.variables[0], func.expr
        return expr.is_rational_function(var)

    @classmethod
    def _rational_case(cls, poly, func):
        """Handle the rational function case. """
        roots = symbols('r:%d' % poly.degree())
        var, expr = func.variables[0], func.expr

        f = sum(expr.subs(var, r) for r in roots)
        p, q = together(f).as_numer_denom()

        domain = QQ[roots]

        p = p.expand()
        q = q.expand()

        try:
            p = Poly(p, domain=domain, expand=False)
        except GeneratorsNeeded:
            p, p_coeff = None, (p,)
        else:
            p_monom, p_coeff = zip(*p.terms())

        try:
            q = Poly(q, domain=domain, expand=False)
        except GeneratorsNeeded:
            q, q_coeff = None, (q,)
        else:
            q_monom, q_coeff = zip(*q.terms())

        coeffs, mapping = symmetrize(p_coeff + q_coeff, formal=True)
        formulas, values = viete(poly, roots), []

        for (sym, _), (_, val) in zip(mapping, formulas):
            values.append((sym, val))

        for i, (coeff, _) in enumerate(coeffs):
            coeffs[i] = coeff.subs(values)

        n = len(p_coeff)

        p_coeff = coeffs[:n]
        q_coeff = coeffs[n:]

        if p is not None:
            p = Poly(dict(zip(p_monom, p_coeff)), *p.gens).as_expr()
        else:
            (p,) = p_coeff

        if q is not None:
            q = Poly(dict(zip(q_monom, q_coeff)), *q.gens).as_expr()
        else:
            (q,) = q_coeff

        return factor(p/q)

    def _hashable_content(self):
        return (self.poly, self.fun)

    @property
    def expr(self):
        return self.poly.as_expr()

    @property
    def args(self):
        return (self.expr, self.fun, self.poly.gen)

    @property
    def free_symbols(self):
        return self.poly.free_symbols | self.fun.free_symbols

    @property
    def is_commutative(self):
        return True

    def doit(self, **hints):
        if not hints.get('roots', True):
            return self

        _roots = roots(self.poly, multiple=True)

        if len(_roots) < self.poly.degree():
            return self
        else:
            return Add(*[ self.fun(r) for r in _roots ])

    def _eval_evalf(self, prec):
        try:
            _roots = self.poly.nroots(n=prec_to_dps(prec))
        except (DomainError, PolynomialError):
            return self
        else:
            return Add(*[ self.fun(r) for r in _roots ])

    def _eval_derivative(self, x):
        var, expr = self.fun.args
        func = Lambda(var, expr.diff(x))
        return self.new(self.poly, func, self.auto)

def bisect(f, a, b, tol):
    """
    Implements bisection. This function is used in RootOf.eval_rational() and
    it needs to be robust.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.polys.rootoftools import bisect
    >>> bisect(lambda x: x**2-1, -10, 0, S(1)/10**2)
    -1025/1024
    >>> bisect(lambda x: x**2-1, -10, 0, S(1)/10**4)
    -131075/131072

    """
    a = sympify(a)
    b = sympify(b)
    fa = f(a)
    fb = f(b)
    if fa * fb >= 0:
        raise ValueError("bisect: f(a) and f(b) must have opposite signs")
    while (b-a > tol):
        c = (a+b)/2
        fc = f(c)
        if (fc == 0): return c # We need to make sure f(c) is not zero below
        if (fa * fc < 0):
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return (a+b)/2
