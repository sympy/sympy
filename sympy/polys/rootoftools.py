"""Implementation of RootOf class and related tools. """

from sympy import Expr, Integer, Real, I, Add, Lambda

from sympy.polys.polytools import Poly

from sympy.polys.rootisolation import (
    dup_isolate_complex_roots_sqf,
    dup_isolate_real_roots_sqf,
)

from sympy.polys.polyroots import (
    roots_linear, roots_quadratic, roots_binomial,
)

from sympy.polys.polyerrors import (
    PolynomialError, DomainError,
)

from sympy.mpmath import (
    mp, mpf, mpc, mpi, findroot,
)

from sympy.utilities import lambdify, any

def dup_minpoly_add(f, g, K):
    """ """
    F = dmp_raise(f, 1, 0, K)
    G = dmp_raise(g, 1, 0, K)

    H = [[-K.one], [K.one, K.zero]]
    F = dmp_compose(F, H, 1, K)

    return dmp_resultant(F, G, 1, K)

def dup_minpoly_sub(f, g, K):
    """ """
    F = dmp_raise(f, 1, 0, K)
    G = dmp_raise(g, 1, 0, K)

    H = [[K.one], [K.one, K.zero]]
    F = dmp_compose(F, H, 1, K)

    return dmp_resultant(F, G, 1, K)

def dup_minpoly_mul(f, g, K):
    """ """
    f, F = reversed(f), []

    for i, c in enumerate(f):
        if not c:
            F.append([])
        else:
            F.append(dup_lshift([c], i, K))

    F = dmp_strip(F)
    G = dmp_raise(g, 1, 0, K)

    return dmp_resultant(F, G, 1, K)

def dup_minpoly_div(f, g, K):
    """ """
    F = dmp_raise(f, 1, 0, K)
    G = dmp_raise(g, 1, 0, K)

    H = [[K.one, K.zero], []]
    F = dmp_compose(F, H, 1, K)

    return dmp_resultant(F, G, 1, K)

def dup_minpoly_pow(f, p, q, K):
    """ """
    d = {(p, 0): -K.one, (0, q): K.one}

    F = dmp_raise(f, 1, 0, K)
    G = dmp_from_dict(d, 1, K)

    return dmp_resultant(F, G, 1, K)

_rootof_trivial_cache = {}

def roots_trivial(poly, radicals=True):
    """Compute roots in linear, quadratic and binomial cases. """
    if poly.degree() == 1:
        return roots_linear(poly)
    else:
        if not radicals:
            return None

        if poly in _rootof_trivial_cache:
            roots = _rootof_trivial_cache[poly]
        else:
            if radicals and poly.degree() == 2:
                roots = roots_quadratic(poly)
            elif radicals and poly.length() == 2 and poly.TC():
                roots = roots_binomial(poly)
            else:
                return None

            _rootof_trivial_cache[poly] = roots

        return roots

_rootof_reals_cache = {}
_rootof_complexes_cache = {}

def _rootof_get_reals_sqf(factor):
    """Compute real isolating intervals for a square-free polynomial. """
    if factor in _rootof_reals_cache:
        real_part = _rootof_reals_cache[factor]
    else:
        _rootof_reals_cache[factor] = real_part = \
            dup_isolate_real_roots_sqf(factor.rep.rep, factor.rep.dom, blackbox=True)

    return real_part

def _rootof_get_complexes_sqf(factor):
    """Compute complex isolating intervals for a square-free polynomial. """
    if factor in _rootof_complexes_cache:
        complex_part = _rootof_complexes_cache[factor]
    else:
        _rootof_complexes_cache[factor] = complex_part = \
            dup_isolate_complex_roots_sqf(factor.rep.rep, factor.rep.dom, blackbox=True)

    return complex_part

def _rootof_get_reals(factors):
    """Compute real isolating intervals for a list of factors. """
    reals = []

    for factor, k in factors:
        real_part = _rootof_get_reals_sqf(factor)
        reals.extend([ (root, factor, k) for root in real_part ])

    return reals

def _rootof_get_complexes(factors):
    """Compute complex isolating intervals for a list of factors. """
    complexes = []

    for factor, k in factors:
        complex_part = _rootof_get_complexes_sqf(factor)
        complexes.extend([ (root, factor, k) for root in complex_part ])

    return complexes

def _rootof_reals_sorted(reals):
    """Make real isolating intervals disjoint and sort roots. """
    cache = {}

    for i, (u, f, k) in enumerate(reals):
        for j, (v, g, m) in enumerate(reals[i+1:]):
            u, v = u.refine_disjoint(v)
            reals[i+j+1] = (v, g, m)

        reals[i] = (u, f, k)

    reals = sorted(reals, key=lambda r: (r[0].a, r[0].b))

    for root, factor, _ in reals:
        if factor in cache:
            cache[factor].append(root)
        else:
            cache[factor] = [root]

    for factor, roots in cache.iteritems():
        _rootof_reals_cache[factor] = roots

    return reals

def _rootof_complexes_sorted(complexes):
    """Make complex isolating intervals disjoint and sort roots. """
    cache = {}

    for i, (u, f, k) in enumerate(complexes):
        for j, (v, g, m) in enumerate(complexes[i+1:]):
            u, v = u.refine_disjoint(v)
            complexes[i+j+1] = (v, g, m)

        complexes[i] = (u, f, k)

    complexes = sorted(complexes, key=lambda r: (r[0].ax, r[0].ay))

    for root, factor, _ in complexes:
        if factor in cache:
            cache[factor].append(root)
        else:
            cache[factor] = [root]

    for factor, roots in cache.iteritems():
        _rootof_complexes_cache[factor] = roots

    return complexes

def _rootof_reals_index(reals, index):
    """Transform ``RootOf`` index concerning real roots. """
    i = 0

    for j, (_, factor, k) in enumerate(reals):
        if index < i + k:
            poly, index = factor, 0

            for _, factor, _ in reals[:j]:
                if factor == poly:
                    index += 1

            return poly, index, None, None
        else:
            i += k

def _rootof_complexes_index(complexes, index):
    """Transform ``RootOf`` index concerning complex roots. """
    index, conjugate, i = index, False, 0

    for j, (_, factor, k) in enumerate(complexes):
        if index < i + 2*k:
            if index >= i + k:
                conjugate = True

            poly, pointer = factor, 0

            for _, factor, _ in complexes[:j]:
                if factor == poly:
                    pointer += 1

            index = len(_rootof_reals_cache[poly])

            if not conjugate:
                index += 2*pointer
            else:
                index += 2*pointer + 1

            return poly, index, pointer, conjugate
        else:
            i += 2*k

def _rootof_data(poly, indices):
    """Construct ``RootOf`` data from a polynomial and indices. """
    (_, factors) = poly.factor_list()

    reals = _rootof_get_reals(factors)
    real_count = sum([ k for _, _, k in reals ])

    if indices is None:
        reals = _rootof_reals_sorted(reals)

        for index in xrange(0, real_count):
            yield _rootof_reals_index(reals, index)
    else:
        if any(index < real_count for index in indices):
            reals = _rootof_reals_sorted(reals)

            for index in indices:
                if index < real_count:
                    yield _rootof_reals_index(reals, index)

        if any(index >= real_count for index in indices):
            complexes = _rootof_get_complexes(factors)
            complexes = _rootof_complexes_sorted(complexes)

            for index in indices:
                if index >= real_count:
                    yield _rootof_complexes_index(complexes, index-real_count)

class RootOf(Expr):
    """Represents ``k``-th root of a univariate polynomial. """

    __slots__ = ['poly', 'index', 'pointer', 'conjugate']

    def __new__(cls, f, indices=None, radicals=True, expand=True):
        """Construct a new ``RootOf`` object for ``k``-th root of ``f``. """
        poly = Poly(f, greedy=False, expand=expand)

        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are supported")

        if poly.degree() <= 0:
            raise PolynomialError("can't construct RootOf object for %s" % f)

        if indices is not None:
            if hasattr(indices, '__iter__'):
                indices, iterable = list(indices), True
            else:
                indices, iterable = [indices], False

            deg = poly.degree()

            for i, index in enumerate(indices):
                if index < -deg or index >= deg:
                    raise IndexError("root index out of [%d, %d] range, got %d" % (-deg, deg-1, index))
                elif index < 0:
                    indices[i] += deg
        else:
            iterable = True

        if not poly.get_domain().is_Exact:
            poly = poly.to_exact()

        roots = roots_trivial(poly, radicals)

        if roots is not None:
            if indices is not None:
                result = [ roots[index] for index in indices ]
            else:
                result = [ root for root in roots if root.is_real ]
        else:
            dom = poly.get_domain()

            if dom.is_QQ:
                _, poly = poly.clear_denoms(convert=True)
            elif not dom.is_ZZ:
                raise DomainError("RootOf is not supported over %s" % dom)

            result = []

            for data in _rootof_data(poly, indices):
                poly, index, pointer, conjugate = data

                roots = roots_trivial(poly, radicals)

                if roots is not None:
                    result.append(roots[index])
                else:
                    result.append(cls._inner_new(poly, index, pointer, conjugate))

        if not iterable:
            return result[0]
        else:
            return result

    @classmethod
    def _inner_new(cls, poly, index, pointer=None, conjugate=None):
        """Construct new ``RootOf`` instance from valid ``RootOf`` data. """
        obj = Expr.__new__(cls)

        obj.poly = poly
        obj.index = index

        if pointer is None:
            obj.pointer = index
        else:
            obj.pointer = pointer

        obj.conjugate = conjugate

        return obj

    def _hashable_content(self):
        return (self.expr, self.index)

    @property
    def expr(self):
        return self.poly.as_basic()

    @property
    def args(self):
        return [self.expr, Integer(self.index)]

    @property
    def is_real(self):
        """Return ``True`` if the root in consideration is real. """
        return self.conjugate is None

    @property
    def is_complex(self):
        """Return ``True`` if the root in consideration is complex. """
        return self.conjugate is not None

    @property
    def is_conjugate(self):
        """Return ``True`` if the root is located in the lower half-plane. """
        return self.is_complex and self.conjugate

    def _get_interval(self):
        """Internal function for retrieving isolation interval from cache. """
        if self.is_real:
            return _rootof_reals_cache[self.poly][self.pointer]
        else:
            return _rootof_complexes_cache[self.poly][self.pointer]

    def _set_interval(self, interval):
        """Internal function for updating isolation interval in cache. """
        if self.is_real:
            _rootof_reals_cache[self.poly][self.pointer] = interval
        else:
            _rootof_complexes_cache[self.poly][self.pointer] = interval

    def _eval_evalf(self, prec):
        """Evaluate this complex root to the given precision. """
        _prec, mp.prec = mp.prec, prec

        try:
            func = lambdify(self.poly.gen, self.expr)
            interval, refined = self._get_interval(), False

            while True:
                if self.is_real:
                    x0 = mpf(str(interval.center))
                else:
                    re, im = interval.center

                    re = mpf(str(re))
                    im = mpf(str(im))

                    x0 = mpc(re, im)

                try:
                    root = findroot(func, x0)
                except ValueError:
                    interval = interval.refine()
                    refined = True
                    continue
                else:
                    if refined:
                        self._set_interval(interval)

                    if self.is_conjugate:
                        root = root.conjugate()

                    break
        finally:
            mp.prec = _prec

        return Real._new(root.real._mpf_, prec) + I*Real._new(root.imag._mpf_, prec)

class RootSum(Expr):
    """Represents a sum of all roots of a univariate polynomial. """

    __slots__ = ['poly', 'func', 'roots']

    def __new__(cls, poly, func=None, formal=True):
        """Construct new ``RootSum`` instance carrying all formal roots of ``poly``. """
        poly = Poly(poly, greedy=False)

        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are supported")

        if func is None:
            func = Lambda(poly.gen, poly.gen)
        elif not hasattr(func, '__call__'):
            raise TypeError("%s is not a callable object" % func)

        (_, factors), terms = poly.factor_list(), []

        for poly, k in factors:
            roots = []

            for i in xrange(0, poly.degree()):
                root = RootOf(poly, i)

                if formal and root.has(RootOf):
                    roots.append(root)
                else:
                    terms.append(k*func(root))

            if formal and roots:
                obj = Expr.__new__(cls)

                obj.poly = poly
                obj.func = func
                obj.roots = roots

                terms.append(k*obj)

        return Add(*terms)

    def _hashable_content(self):
        return (self.expr, self.func)

    @property
    def expr(self):
        return self.poly.as_basic()

    @property
    def args(self):
        return [self.expr, self.func]

    def doit(self, **hints):
        if hints.get('roots', True):
            return Add(*map(self.func, self.roots))
        else:
            return self

