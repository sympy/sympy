"""Exact real algebraic sample points for cylindrical algebraic decomposition.

The lifting phase of CAD needs to evaluate polynomials at points whose
coordinates are real algebraic numbers, to find the real roots of the
resulting univariate polynomials and to compare them exactly. Here a sample
point is represented by a single real algebraic number `\\theta`, given as a
:class:`~.ComplexRootOf` with its isolating interval, together with the
coordinates as elements of the field `\\mathbb{Q}(\\theta)`. Signs and
comparisons are decided by interval arithmetic, refining the isolating
interval of `\\theta` as needed, and the field is extended by a primitive
element whenever a new algebraic coordinate is added.
"""
from __future__ import annotations

from itertools import count

from sympy.core.expr import Expr
from sympy.core.numbers import Rational
from sympy.core.symbol import Dummy
from sympy.polys.densetools import dup_eval
from sympy.polys.domains import QQ
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.polytools import Poly
from sympy.polys.rootoftools import ComplexRootOf


def _split(r):
    """Split a real algebraic number into a rational factor and a
    :class:`~.ComplexRootOf` (or ``None`` for a rational number).

    Real roots are represented by :class:`~.Rational` numbers, by real
    :class:`~.ComplexRootOf` roots or by products of a rational number and
    such a root, which is the form :func:`~.real_roots` may return.
    """
    if isinstance(r, ComplexRootOf):
        return QQ.one, r
    if isinstance(r, Expr) and r.is_Mul and len(r.args) == 2 and \
            r.args[0].is_Rational and isinstance(r.args[1], ComplexRootOf):
        return QQ.convert(r.args[0]), r.args[1]
    return QQ.convert(r), None


def _bounds(r):
    """Rational lower and upper bounds of the real algebraic number ``r``."""
    c, root = _split(r)
    if root is None:
        return c, c
    interval = root._get_interval()
    a, b = c*QQ.convert(interval.a), c*QQ.convert(interval.b)
    return (a, b) if a <= b else (b, a)


def _refine(r):
    """Refine the isolating interval of the real algebraic number ``r``."""
    _, root = _split(r)
    if root is not None:
        root._set_interval(root._get_interval().refine())


def _minpoly(r):
    """Minimal polynomial over ZZ of the real algebraic number ``r``, which
    must not be rational."""
    c, root = _split(r)
    p = root.poly
    if c != 1:
        x = p.gen
        d = p.degree()
        p = Poly(p.as_expr().subs(x, x/c.numerator*c.denominator)*c.numerator**d, x)
        p = p.primitive()[1]
    return p


def compare_real(a, b):
    """Compare two real algebraic numbers exactly.

    ``a`` and ``b`` must be :class:`~.Rational` numbers, real
    :class:`~.ComplexRootOf` roots or products of the two. Returns ``-1``,
    ``0`` or ``1``.

    Examples
    ========

    >>> from sympy import CRootOf, Rational
    >>> from sympy.abc import x
    >>> from sympy.polys.cad.samplepoints import compare_real
    >>> compare_real(CRootOf(x**2 - 2, 1), Rational(3, 2))
    -1
    >>> compare_real(CRootOf(x**2 - 2, 1), CRootOf(x**3 - 2, 0))
    1
    """
    if a == b:
        return 0
    while True:
        alo, ahi = _bounds(a)
        blo, bhi = _bounds(b)
        if ahi < blo:
            return -1
        if bhi < alo:
            return 1
        _refine(a)
        _refine(b)


def _floor(q):
    return q.numerator // q.denominator


def _floor_scaled(a, k):
    """Exact value of `\\lfloor 2^k a \\rfloor` for a real algebraic ``a``."""
    scale = QQ(2)**k
    while True:
        lo, hi = _bounds(a)
        flo, fhi = _floor(lo*scale), _floor(hi*scale)
        if flo == fhi:
            return flo
        _refine(a)


def simplest_between(lo, hi, lo_open=True, hi_open=True):
    """The rational number with the smallest denominator in the interval
    with rational endpoints ``lo`` and ``hi``. The endpoints are excluded
    unless ``lo_open`` or ``hi_open`` is ``False``; ``hi`` may be ``None``
    for `+\\infty`.

    Examples
    ========

    >>> from sympy import QQ
    >>> from sympy.polys.cad.samplepoints import simplest_between
    >>> simplest_between(QQ(1, 3), QQ(1, 2))
    2/5
    >>> simplest_between(QQ(1, 3), QQ(1, 2), hi_open=False)
    1/2
    >>> simplest_between(QQ(5, 2), QQ(3))
    8/3
    """
    if hi is None:
        if lo < 0 or (lo == 0 and not lo_open):
            return QQ.zero
        n = _floor(lo)
        return QQ(n if n == lo and not lo_open else n + 1)
    if lo > hi or (lo == hi and (lo_open or hi_open)):
        raise ValueError("empty interval")
    if (lo < 0 or (lo == 0 and not lo_open)) and (hi > 0 or (hi == 0 and not hi_open)):
        return QQ.zero
    if hi < 0 or (hi == 0 and hi_open):
        return -simplest_between(-hi, -lo, hi_open, lo_open)
    # now 0 <= lo < hi, with lo excluded if it is 0
    n = _floor(lo)
    c = n if (n == lo and not lo_open) else n + 1
    if c < hi or (c == hi and not hi_open):
        return QQ(c)
    # no integer is admissible: both bounds lie in [n, n + 1] and the
    # result is n + 1/q with q in the reciprocal interval
    inner_hi = 1/(lo - n) if lo > n else None
    inner = simplest_between(1/(hi - n), inner_hi, hi_open, lo_open)
    return n + 1/inner


def _dyadic_bounds(r, k):
    """Rational bounds of the real algebraic number ``r`` at the scale
    `2^{-k}`, and whether they are attained.

    For a rational number both bounds are the number itself and are
    excluded; for an irrational one they are the dyadic numbers just below
    and above it, and are included.
    """
    c, root = _split(r)
    if root is None:
        return c, c, True, True
    scale = QQ(2)**k
    lo = QQ(_floor_scaled(r, k))/scale
    return lo, lo + 1/scale, False, False


def rational_between(a, b):
    """A simple rational number strictly between the real algebraic
    numbers ``a < b``.

    The result only depends on ``a`` and ``b``: it is the rational number
    with the smallest denominator between the dyadic approximations of the
    irrational endpoints, at the coarsest scale that separates them.

    Examples
    ========

    >>> from sympy import CRootOf
    >>> from sympy.abc import x
    >>> from sympy.polys.cad.samplepoints import rational_between
    >>> rational_between(CRootOf(x**2 - 2, 0), CRootOf(x**2 - 2, 1))
    0
    >>> rational_between(1, CRootOf(x**2 - 2, 1))
    5/4
    >>> rational_between(-1, 1)
    0
    """
    if compare_real(a, b) >= 0:
        raise ValueError("%s is not smaller than %s" % (a, b))
    for k in count(0):
        _, hi, _, hi_open = _dyadic_bounds(a, k)
        lo, _, lo_open, _ = _dyadic_bounds(b, k)
        if hi < lo or (hi == lo and not (hi_open or lo_open)):
            return Rational(simplest_between(hi, lo, hi_open, lo_open))


def rational_below(a):
    """A simple rational number smaller than the real algebraic number ``a``.

    >>> from sympy import CRootOf
    >>> from sympy.abc import x
    >>> from sympy.polys.cad.samplepoints import rational_below
    >>> rational_below(CRootOf(x**2 - 2, 1))
    0
    >>> rational_below(-3)
    -4
    """
    if compare_real(a, Rational(0)) > 0:
        return Rational(0)
    fl = _floor_scaled(a, 0)
    if fl == a:
        fl -= 1
    return Rational(fl)


def rational_above(a):
    """A simple rational number greater than the real algebraic number ``a``.

    >>> from sympy import CRootOf
    >>> from sympy.abc import x
    >>> from sympy.polys.cad.samplepoints import rational_above
    >>> rational_above(CRootOf(x**2 - 2, 1))
    2
    """
    if compare_real(a, Rational(0)) < 0:
        return Rational(0)
    return Rational(_floor_scaled(a, 0) + 1)


def _sign_in_field(theta, a):
    """Sign of the element ``a`` of `\\mathbb{Q}(\\theta)`, with
    ``theta`` a real :class:`~.ComplexRootOf`."""
    if not a:
        return 0
    coeffs = a.to_list()
    while True:
        lo, hi = _bounds(theta)
        vlo = vhi = QQ.zero
        for c in coeffs:
            products = [vlo*lo, vlo*hi, vhi*lo, vhi*hi]
            vlo, vhi = min(products) + c, max(products) + c
        if vlo > 0:
            return 1
        if vhi < 0:
            return -1
        _refine(theta)


def _horner(coeffs, a, K):
    """Evaluate the polynomial with rational ``coeffs`` at the element
    ``a`` of the field ``K``."""
    result = K.zero
    for c in coeffs:
        result = result*a + K.convert(c, QQ)
    return result


def _join(theta, beta):
    """Primitive element of `\\mathbb{Q}(\\theta, \\beta)`.

    Returns ``(gamma, K, theta_K, beta_K)`` where ``gamma`` is a real
    :class:`~.ComplexRootOf` generating the field, ``K`` is the algebraic
    field ``QQ<gamma>`` and ``theta_K``, ``beta_K`` are ``theta`` and
    ``beta`` as elements of ``K``.
    """
    z, t = Dummy('z'), Dummy('t')
    m1, m2 = _minpoly(theta), _minpoly(beta)
    for k in count(1):
        # the roots of N are the sums of a conjugate of theta and k times a
        # conjugate of beta: if they are all distinct then gamma is
        # primitive
        m1z = Poly(m1.as_expr().subs(m1.gen, z - k*t), t, z)
        m2t = Poly(m2.as_expr().subs(m2.gen, t), t, z)
        N = m2t.resultant(m1z)
        if N.is_sqf:
            break
    candidates = []
    for factor, _ in N.factor_list()[1]:
        candidates.extend(factor.real_roots(radicals=False))
    while True:
        tlo, thi = _bounds(theta)
        blo, bhi = _bounds(beta)
        lo, hi = tlo + k*blo, thi + k*bhi
        alive = []
        for r in candidates:
            rlo, rhi = _bounds(r)
            if not (rhi < lo or hi < rlo):
                alive.append(r)
        if len(alive) == 1:
            break
        # the candidates come from different factors, so their isolating
        # intervals may overlap: refine them as well
        _refine(theta)
        _refine(beta)
        for r in alive:
            _refine(r)
    gamma = alive[0]
    if _split(gamma)[1] is None:
        raise PolynomialError("unexpected rational primitive element")
    K = QQ.algebraic_field(gamma)
    # beta is the only common root of m2(t) and m1(gamma - k*t)
    u = Poly.from_list([K.convert(-k), K.unit], t, domain=K)
    f = Poly(0, t, domain=K)
    for c in m1.all_coeffs():
        f = f*u + c
    g = Poly(m2.as_expr().subs(m2.gen, t), t, domain=K).gcd(f)
    if g.degree() != 1:
        raise PolynomialError("failed to compute a primitive element")
    g1, g0 = g.rep.to_list()
    beta_K = -g0/g1
    theta_K = K.unit - k*beta_K
    return gamma, K, theta_K, beta_K


class SamplePoint:
    """A point with real algebraic coordinates.

    The coordinates are elements of the field ``QQ<theta>`` (or ``QQ``) where
    ``theta`` is a real :class:`~.ComplexRootOf`. Points are built by
    extending the origin one coordinate at a time with :meth:`extend`.

    Examples
    ========

    >>> from sympy import CRootOf, Rational
    >>> from sympy.abc import x, y
    >>> from sympy.polys.cad.samplepoints import SamplePoint
    >>> q = SamplePoint().extend(CRootOf(x**2 - 2, 1))
    >>> q.real_roots(x**2 + y**2 - 1, [x, y])
    []
    >>> q.real_roots(x**2 - y**2 - 1, [x, y])
    [-1, 1]
    >>> p = q.extend(Rational(1, 2))
    >>> p
    SamplePoint(CRootOf(x**2 - 2, 1), 1/2)
    >>> p.sign(x**2 + y**2 - 2, [x, y])
    1
    >>> p.sign(2*y - 1, [x, y])
    0
    """

    __slots__ = ('field', 'theta', 'coords')

    def __init__(self, field=QQ, theta=None, coords=()):
        self.field = field
        self.theta = theta
        self.coords = tuple(coords)

    def __len__(self):
        return len(self.coords)

    def __repr__(self):
        return "SamplePoint(%s)" % ", ".join(str(c) for c in self.as_exprs())

    def as_exprs(self):
        """The coordinates as SymPy expressions."""
        return tuple(self.field.to_sympy(c) for c in self.coords)

    @property
    def degree(self):
        """Degree of the field of the coordinates over the rationals."""
        if self.theta is None:
            return 1
        return _minpoly(self.theta).degree()

    def extend(self, root):
        """The point with the coordinate ``root`` appended.

        ``root`` is a :class:`~.Rational`, a real :class:`~.ComplexRootOf`
        or a product of the two.
        """
        if _split(root)[1] is not None:
            if root == self.theta:
                return SamplePoint(self.field, self.theta,
                                   self.coords + (self.field.unit,))
            if self.theta is None:
                K = QQ.algebraic_field(root)
                coords = [K.convert(c, QQ) for c in self.coords]
                return SamplePoint(K, root, coords + [K.unit])
            gamma, K, theta_K, beta_K = _join(self.theta, root)
            coords = [_horner(c.to_list(), theta_K, K) for c in self.coords]
            return SamplePoint(K, gamma, coords + [beta_K])
        value = self.field.convert(QQ.from_sympy(Rational(root)), QQ)
        return SamplePoint(self.field, self.theta, self.coords + (value,))

    def _evaluate(self, poly, gens):
        """``poly`` in the generators ``gens`` with the leading ones
        replaced by the coordinates: a polynomial in the remaining generators
        over the field, or an element of the field if none remains."""
        gens = tuple(gens)
        n = len(self.coords)
        if len(gens) < n:
            raise ValueError("the point has more coordinates than generators")
        if not gens:
            return QQ.convert(poly.as_expr() if isinstance(poly, Poly) else poly)
        if not isinstance(poly, Poly) or poly.gens != gens:
            poly = Poly(poly, *gens)
        f = poly.set_domain(self.field)
        for gen, value in zip(gens[:n - 1], self.coords):
            f = f.eval(gen, value)
        if n:
            last = self.coords[n - 1]
            if len(gens) == n:
                return dup_eval(f.rep.to_list(), last, self.field)
            f = f.eval(gens[n - 1], last)
        return f

    def sign(self, poly, gens):
        """Exact sign of the polynomial ``poly`` in the generators ``gens``
        at the point. There must be as many generators as coordinates."""
        if len(gens) != len(self.coords):
            raise ValueError("expected %d generators" % len(self.coords))
        value = self._evaluate(poly, gens)
        if self.theta is None:
            return (value > 0) - (value < 0)
        return _sign_in_field(self.theta, value)

    def real_roots(self, poly, gens):
        """Sorted distinct real roots of ``poly`` in the last generator of
        ``gens``, after substituting the point for the other generators.

        Returns ``None`` if the polynomial vanishes identically at the
        point.
        """
        if len(gens) != len(self.coords) + 1:
            raise ValueError("expected %d generators" % (len(self.coords) + 1))
        f = self._evaluate(poly, gens)
        if f.is_zero:
            return None
        if f.degree() <= 0:
            return []
        roots = [r for r, _ in f.real_roots(multiple=False, radicals=False)]
        return sorted(roots, key=_SortKey)


class _SortKey:
    __slots__ = ('root',)

    def __init__(self, root):
        self.root = root

    def __lt__(self, other):
        return compare_real(self.root, other.root) < 0
