"""Implementation of RootOf class and related tools. """

from sympy import Basic, Integer, Real, I

from sympy.polys.polytools import Poly

from sympy.polys.rootisolation import (
    RealInterval, ComplexInterval,
    dup_isolate_complex_roots_sqf,
    dup_isolate_real_roots_sqf,
)

from sympy.polys.polyerrors import (
    PolynomialError, DomainError,
)

from sympy.mpmath import (
    mp, mpf, mpc, mpi, findroot,
)

from sympy.utilities import lambdify

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

_rootof_reals_cache = {}
_rootof_complexes_cache = {}

class RootOf(Basic):
    """Represents ``k``-th root of a univariate polynomial. """

    __slots__ = ['poly', 'index', 'interval']

    def __new__(cls, f, index, expand=True):
        """Construct a new ``RootOf`` object for ``k``-th root of ``f``. """
        poly = Poly(f, greedy=False, expand=expand)

        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are supported")

        dom = poly.get_domain()

        if dom.is_QQ:
            _, poly = poly.clear_denoms(convert=True)
        elif not dom.is_ZZ:
            raise DomainError("RootOf is not supported over %s" % dom)

        deg = poly.degree()

        if index < -deg or index >= deg:
            raise IndexError("root index out of (%d, %d) range, got %d" % (-deg, deg-1, index))
        elif index < 0:
            index += deg

        poly, index, interval = cls._inner_init(poly, index)

        if poly.degree() == 1:
            return -poly.nth(0)/poly.nth(1)

        obj = Basic.__new__(cls)

        obj.poly = poly
        obj.index = index
        obj.interval = interval

        return obj

    @classmethod
    def _inner_init(cls, poly, index):
        """Decompose the input polynomial and compute ``RootOf`` data. """
        (_, factors), reals = poly.factor_list(), []

        for f, k in factors:
            if f in _rootof_reals_cache:
                real_part = _rootof_reals_cache[f]
            else:
                _rootof_reals_cache[f] = real_part = \
                    dup_isolate_real_roots_sqf(f.rep.rep, f.rep.dom, blackbox=True)

            reals.extend([ (r, f, k) for r in real_part ])

        real_count = sum([ k for _, _, k in reals ])

        if index < real_count:
            for i, (u, f, k) in enumerate(reals):
                for j, (v, g, m) in enumerate(reals[i+1:]):
                    u, v = u.refine_disjoint(v)
                    reals[i+j+1] = (v, g, m)

                reals[i] = (u, f, k)

            reals = sorted(reals, key=lambda r: (r[0].a, r[0].b))

            cache, poly, i = {}, None, 0

            for j, (root, factor, k) in enumerate(reals):
                if poly is None:
                    if index < i + k:
                        poly, index, interval = factor, 0, root

                        for _, f, _ in reals[:j]:
                            if f == factor:
                                index += 1
                    else:
                        i += k

                if factor in cache:
                    cache[factor].append(root)
                else:
                    cache[factor] = [root]

            for factor, reals in cache.iteritems():
                _rootof_reals_cache[factor] = reals
        else:
            index, complexes = index - real_count, []

            for f, k in factors:
                if f in _rootof_complexes_cache:
                    complex_part = _rootof_complexes_cache[f]
                else:
                    _rootof_complexes_cache[f] = complex_part = \
                        dup_isolate_complex_roots_sqf(f.rep.rep, f.rep.dom, blackbox=True)

                complexes.extend([ (r, f, k) for r in complex_part ])

            for i, (u, f, k) in enumerate(complexes):
                for j, (v, g, m) in enumerate(complexes[i+1:]):
                    u, v = u.refine_disjoint(v)
                    complexes[i+j+1] = (v, g, m)

                complexes[i] = (u, f, k)

            complexes = sorted(complexes, key=lambda r: (r[0].ax, r[0].ay))

            cache, poly, conjugate, i = {}, None, False, 0

            for j, (root, factor, k) in enumerate(complexes):
                if poly is None:
                    if index < i + 2*k:
                        if index >= i + k:
                            conjugate = True

                        poly, index = factor, len(_rootof_reals_cache[factor])

                        if conjugate:
                            index += 1

                        for _, f, _ in complexes[:j]:
                            if f == factor:
                                index += 2

                        if not conjugate:
                            interval = root
                        else:
                            interval = root.conjugate()
                    else:
                        i += 2*k

                if factor in cache:
                    cache[factor].append(root)
                else:
                    cache[factor] = [root]

            for factor, complexes in cache.iteritems():
                _rootof_complexes_cache[factor] = complexes

        return poly, index, interval

    @property
    def expr(self):
        return self.poly.as_basic()

    @property
    def args(self):
        return [self.expr, Integer(self.index)]

    @property
    def is_real(self):
        """Return ``True`` if the root in consideration is real. """
        return isinstance(self.interval, RealInterval)

    @property
    def is_complex(self):
        """Return ``True`` if the root in consideration is complex. """
        return isinstance(self.interval, ComplexInterval)

    @property
    def is_conjugate(self):
        """Return ``True`` if the root is located in the lower half-plane. """
        return self.is_complex and self.interval.conj

    def _eval_evalf(self, prec):
        """Evaluate this complex root to the given precision. """
        _prec, mp.prec = mp.prec, prec

        try:
            func = lambdify(self.poly.gen, self.expr)
            interval, refined = self.interval, False

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
                    break
        finally:
            mp.prec = _prec

        return Real._new(root.real._mpf_, prec) + I*Real._new(root.imag._mpf_, prec)
