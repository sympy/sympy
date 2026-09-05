from __future__ import annotations

from sympy.core.numbers import Rational
from sympy.functions.elementary.integers import floor
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys.domains import QQ
from sympy.polys.polytools import Poly
from sympy.polys.rootoftools import CRootOf
from sympy.polys.cad.samplepoints import (SamplePoint, compare_real,
    simplest_between, rational_between, rational_below, rational_above,
    _floor_scaled, _join, _sign_in_field, _minpoly)
from sympy.testing.pytest import raises
from sympy.abc import x, y, z, w


def _roots(f):
    return Poly(f, x).real_roots(radicals=False)


def test_simplest_between():
    assert simplest_between(QQ(1, 3), QQ(1, 2)) == QQ(2, 5)
    assert simplest_between(QQ(5, 2), QQ(3)) == QQ(8, 3)
    assert simplest_between(QQ(-1), QQ(1)) == 0
    assert simplest_between(QQ(2), QQ(5)) == 3
    assert simplest_between(QQ(0), QQ(1)) == QQ(1, 2)
    assert simplest_between(QQ(1), QQ(2)) == QQ(3, 2)
    assert simplest_between(QQ(-3, 2), QQ(-1)) == QQ(-4, 3)
    assert simplest_between(QQ(7, 3), None) == 3
    assert simplest_between(QQ(-7, 3), None) == 0
    assert simplest_between(QQ(3), None) == 4
    assert simplest_between(QQ(99, 100), QQ(101, 100)) == 1
    assert simplest_between(QQ(1, 3), QQ(1, 2), hi_open=False) == QQ(1, 2)
    assert simplest_between(QQ(1, 3), QQ(1, 2), lo_open=False) == QQ(1, 3)
    assert simplest_between(QQ(0), QQ(1), lo_open=False) == 0
    assert simplest_between(QQ(-2), QQ(-1), hi_open=False) == -1
    assert simplest_between(QQ(-2), QQ(-1), lo_open=False) == -2
    assert simplest_between(QQ(3), None, lo_open=False) == 3
    assert simplest_between(QQ(1, 2), QQ(1, 2), False, False) == QQ(1, 2)
    assert simplest_between(QQ(7, 5), QQ(3, 2), hi_open=False) == QQ(3, 2)
    raises(ValueError, lambda: simplest_between(QQ(1), QQ(1), lo_open=False))
    assert simplest_between(QQ(100, 101), QQ(101, 102)) == QQ(2011, 2030) or \
        simplest_between(QQ(100, 101), QQ(101, 102)).denominator <= 2030
    raises(ValueError, lambda: simplest_between(QQ(1), QQ(1)))
    raises(ValueError, lambda: simplest_between(QQ(2), QQ(1)))

    # the result has the smallest denominator among all rationals in the
    # interval, checked by brute force
    for lo, hi in [(QQ(1, 3), QQ(1, 2)), (QQ(5, 2), QQ(3)), (QQ(-3, 2), QQ(-1)),
                   (QQ(7, 10), QQ(8, 10)), (QQ(-5, 7), QQ(-2, 3))]:
        q = simplest_between(lo, hi)
        assert lo < q < hi
        for d in range(1, q.denominator):
            for n in range(int(lo*d) - 1, int(hi*d) + 2):
                assert not (lo < QQ(n, d) < hi)


def test_compare_real():
    s2, ms2 = CRootOf(x**2 - 2, 1), CRootOf(x**2 - 2, 0)
    c2 = CRootOf(x**3 - 2, 0)
    assert compare_real(s2, Rational(3, 2)) == -1
    assert compare_real(Rational(3, 2), s2) == 1
    assert compare_real(s2, c2) == 1
    assert compare_real(c2, s2) == -1
    assert compare_real(s2, s2) == 0
    assert compare_real(ms2, s2) == -1
    assert compare_real(Rational(1), Rational(2)) == -1
    assert compare_real(Rational(2), Rational(2)) == 0
    assert compare_real(Rational(1), s2) == -1
    assert compare_real(Rational(2), s2) == 1
    # close roots of different polynomials
    a = CRootOf(x**2 - 2, 1)
    b = CRootOf(x**4 - 4, 1)  # this is sqrt(2) too but CRootOf canonicalizes
    assert a == b
    c = CRootOf(1000000*x**2 - 2000001, 1)
    assert compare_real(a, c) == -1
    assert compare_real(c, a) == 1
    # roots with a rational factor
    r0, r1 = Poly(2*x**2 - 9, x).real_roots(radicals=False)
    assert r1 == 3*CRootOf(2*x**2 - 1, 1)
    assert compare_real(r0, r1) == -1
    assert compare_real(r1, 2) == 1
    assert compare_real(r1, Rational(3, 2)) == 1
    assert compare_real(r1, CRootOf(x**2 - 5, 1)) == -1
    assert compare_real(r0, -2) == -1
    assert _minpoly(r0).all_coeffs() == _minpoly(r1).all_coeffs() == [2, 0, -9]
    assert _minpoly(r0).gens == r0.args[1].poly.gens
    assert _minpoly(CRootOf(x**3 - 2, 0)).all_coeffs() == [1, 0, 0, -2]
    assert rational_between(r0, r1) == 0
    assert rational_between(2, r1) == Rational(33, 16)
    assert rational_below(r0) == -3
    assert rational_above(r1) == 3


def test_floor_scaled():
    s2 = CRootOf(x**2 - 2, 1)
    for k in range(6):
        assert _floor_scaled(s2, k) == floor(2**k*sqrt(2))
    ms2 = CRootOf(x**2 - 2, 0)
    for k in range(6):
        assert _floor_scaled(ms2, k) == floor(-2**k*sqrt(2))
    assert _floor_scaled(Rational(7, 2), 0) == 3
    assert _floor_scaled(Rational(7, 2), 1) == 7
    assert _floor_scaled(Rational(-7, 2), 0) == -4


def test_rational_between():
    s2, ms2 = CRootOf(x**2 - 2, 1), CRootOf(x**2 - 2, 0)
    assert rational_between(ms2, s2) == 0
    assert rational_between(1, s2) == Rational(5, 4)
    assert rational_between(Rational(1, 3), Rational(1, 2)) == Rational(2, 5)
    assert rational_between(-2, ms2) == Rational(-3, 2)
    assert rational_between(ms2, -1) == Rational(-5, 4)
    assert rational_between(-1, 1) == 0
    assert rational_between(Rational(9, 10), s2) == 1
    assert rational_between(ms2, 1) == 0
    assert rational_between(Rational(-1, 2), Rational(1, 3)) == 0
    assert rational_between(Rational(1, 3), Rational(1, 2)) == Rational(2, 5)
    assert rational_between(3, 4) == Rational(7, 2)
    raises(ValueError, lambda: rational_between(s2, 1))
    raises(ValueError, lambda: rational_between(s2, s2))

    # the result does not depend on the state of the isolating intervals
    s2.evalf(60)
    assert rational_between(1, s2) == Rational(5, 4)
    assert rational_between(ms2, s2) == 0

    roots = sorted(_roots((x**2 - 2)*(x**3 - 2)*(x**2 - 3)*(4*x**2 - 5)*(x - 1)),
                   key=lambda r: r.evalf())
    for a, b in zip(roots, roots[1:]):
        q = rational_between(a, b)
        assert q.is_Rational
        assert compare_real(a, q) == -1 and compare_real(q, b) == -1
        assert q.q <= 16


def test_rational_below_above():
    s2, ms2 = CRootOf(x**2 - 2, 1), CRootOf(x**2 - 2, 0)
    assert rational_below(s2) == 0
    assert rational_above(s2) == 2
    assert rational_below(ms2) == -2
    assert rational_above(ms2) == 0
    assert rational_below(2) == 0
    assert rational_above(2) == 3
    assert rational_below(-3) == -4
    assert rational_above(-3) == 0
    assert rational_below(Rational(-7, 2)) == -4
    assert rational_above(Rational(7, 2)) == 4
    assert rational_below(Rational(1, 2)) == 0
    assert rational_above(Rational(-1, 2)) == 0
    assert rational_below(0) == -1
    assert rational_above(0) == 1
    r = CRootOf(x**3 - 100, 0)
    assert rational_below(r) == 0
    assert rational_above(r) == 5
    r = CRootOf(x**3 + 100, 0)
    assert rational_below(r) == -5
    assert rational_above(r) == 0


def test_sign_in_field():
    s2 = CRootOf(x**2 - 2, 1)
    K = QQ.algebraic_field(s2)
    assert _sign_in_field(s2, K.zero) == 0
    assert _sign_in_field(s2, K.unit) == 1
    assert _sign_in_field(s2, -K.unit) == -1
    # sqrt(2) - 1.4142 is positive although small
    a = K.unit - K.convert(QQ(14142, 10000))
    assert _sign_in_field(s2, a) == 1
    a = K.unit - K.convert(QQ(141421356237, 10**11))
    assert _sign_in_field(s2, a) == 1
    a = K.unit - K.convert(QQ(141421356238, 10**11))
    assert _sign_in_field(s2, a) == -1
    # elements whose leading representation coefficient has the wrong sign
    a = K.unit - K.convert(QQ(2))
    assert _sign_in_field(s2, a) == -1
    a = -K.unit + K.convert(QQ(2))
    assert _sign_in_field(s2, a) == 1


def test_join():
    s2, ms2 = CRootOf(x**2 - 2, 1), CRootOf(x**2 - 2, 0)
    s3 = CRootOf(x**2 - 3, 1)
    g, K, tK, bK = _join(s2, s3)
    assert g.poly.degree() == 4
    assert abs(K.to_sympy(tK).evalf(20) - sqrt(2).evalf(20)) < 1e-15
    assert abs(K.to_sympy(bK).evalf(20) - sqrt(3).evalf(20)) < 1e-15
    assert abs(g.evalf(20) - sqrt(2).evalf(20) - sqrt(3).evalf(20)) < 1e-15
    assert _sign_in_field(g, tK**2 - K.convert(QQ(2))) == 0
    assert _sign_in_field(g, bK**2 - K.convert(QQ(3))) == 0
    # conjugates and elements of the same field give no bigger field
    g, K, tK, bK = _join(s2, ms2)
    assert g.poly.degree() == 2
    assert K.to_sympy(tK) == -ms2 and K.to_sympy(bK) == ms2
    assert _sign_in_field(g, tK + bK) == 0
    g, K, tK, bK = _join(s2, CRootOf(x**2 - 8, 1))
    assert g.poly.degree() == 2
    assert K.to_sympy(tK)*2 == K.to_sympy(bK)
    c2 = CRootOf(x**3 - 2, 0)
    g, K, tK, bK = _join(c2, s2)
    assert g.poly.degree() == 6
    assert abs(K.to_sympy(tK).evalf(20) - c2.evalf(20)) < 1e-15
    assert abs(K.to_sympy(bK).evalf(20) - s2.evalf(20)) < 1e-15


def test_sample_point():
    origin = SamplePoint()
    assert len(origin) == 0 and origin.degree == 1
    assert origin.sign(3, []) == 1
    assert origin.sign(Rational(-1, 2), []) == -1
    assert origin.sign(0, []) == 0
    assert origin.real_roots(x**3 - x, [x]) == [-1, 0, 1]
    assert origin.real_roots(x**2 + 1, [x]) == []
    assert origin.real_roots(Poly(x**2 - 2, x), [x]) == [
        CRootOf(x**2 - 2, 0), CRootOf(x**2 - 2, 1)]
    assert origin.real_roots(Poly(0, x), [x]) is None
    assert origin.real_roots(7, [x]) == []
    raises(ValueError, lambda: origin.sign(x, [x]))
    raises(ValueError, lambda: origin.real_roots(x, [x, y]))

    p = origin.extend(Rational(1, 2))
    assert p.field == QQ and p.theta is None and len(p) == 1
    assert p.as_exprs() == (Rational(1, 2),)
    assert p.sign(2*x - 1, [x]) == 0
    assert p.sign(x - 1, [x]) == -1
    assert p.sign(Poly(x**2, x), [x]) == 1
    assert p.real_roots(y**2 - x, [x, y]) == [-sqrt(2)/2, sqrt(2)/2] or \
        p.real_roots(y**2 - x, [x, y]) == [CRootOf(2*x**2 - 1, 0), CRootOf(2*x**2 - 1, 1)]
    assert p.real_roots(x*y - y, [x, y]) == [0]
    assert p.real_roots((2*x - 1)*y, [x, y]) is None

    s2 = CRootOf(x**2 - 2, 1)
    q = origin.extend(s2)
    assert q.degree == 2 and q.theta == s2
    assert q.as_exprs() == (s2,)
    assert repr(q) == "SamplePoint(CRootOf(x**2 - 2, 1))"
    assert q.sign(x**2 - 2, [x]) == 0
    assert q.sign(x - 1, [x]) == 1
    assert q.sign(x - 2, [x]) == -1
    assert q.real_roots(x**2 + y**2 - 1, [x, y]) == []
    assert q.real_roots(x**2 - y**2 - 1, [x, y]) == [-1, 1]
    assert q.real_roots(x**2 - 2, [x, y]) is None
    assert q.real_roots(x**2 - 1, [x, y]) == []
    r = q.real_roots(y**2 - x, [x, y])
    assert r == [CRootOf(x**4 - 2, 0), CRootOf(x**4 - 2, 1)]

    # the same algebraic number again does not extend the field
    qq = q.extend(s2)
    assert qq.degree == 2 and qq.as_exprs() == (s2, s2)
    assert qq.sign(x - y, [x, y]) == 0
    # a rational coordinate neither
    p2 = q.extend(Rational(1, 2))
    assert p2.degree == 2
    assert repr(p2) == "SamplePoint(CRootOf(x**2 - 2, 1), 1/2)"
    assert p2.sign(x**2 + y**2 - 2, [x, y]) == 1
    assert p2.sign(2*y - 1, [x, y]) == 0
    assert p2.sign(x*y - 1, [x, y]) == -1

    # extension by a root over an algebraic point
    q2 = q.extend(r[1])
    assert q2.degree == 4
    e1, e2 = q2.as_exprs()
    assert abs(e1.evalf(20) - sqrt(2).evalf(20)) < 1e-15
    assert abs(e2.evalf(20) - 2**Rational(1, 4)) < 1e-15
    assert q2.sign(y**2 - x, [x, y]) == 0
    assert q2.sign(y - x, [x, y]) == -1
    assert q2.sign(y**4 - 2, [x, y]) == 0
    r3 = q2.real_roots(z**2 - x - y, [x, y, z])
    assert len(r3) == 2 and compare_real(r3[0], r3[1]) == -1
    q3 = q2.extend(r3[1])
    assert q3.degree == 8
    assert q3.sign(z**2 - x - y, [x, y, z]) == 0
    assert q3.sign(z - 1, [x, y, z]) == 1
    assert q3.sign(z - 2, [x, y, z]) == -1
    assert q3.real_roots(0*w, [x, y, z, w]) is None
    assert q3.real_roots(w**2 + 1, [x, y, z, w]) == []
    assert q3.real_roots(w**2 - z**2, [x, y, z, w]) == [-r3[1], r3[1]] or \
        len(q3.real_roots(w**2 - z**2, [x, y, z, w])) == 2
    raises(ValueError, lambda: q3.sign(x, [x]))
    raises(ValueError, lambda: q3.real_roots(w, [x, y]))

    # distinct roots are sorted and duplicates removed
    roots = q2.real_roots((z**2 - x)*(z - y)*(z**2 - 3), [x, y, z])
    assert roots == [CRootOf(x**2 - 3, 0), CRootOf(x**4 - 2, 0),
                     CRootOf(x**4 - 2, 1), CRootOf(x**2 - 3, 1)]

    # coordinates given as a rational multiple of a root
    r0, r1 = Poly(2*x**2 - 9, x).real_roots(radicals=False)
    p = origin.extend(r1)
    assert p.degree == 2 and p.theta == r1
    assert p.sign(2*x**2 - 9, [x]) == 0
    assert p.sign(x - 2, [x]) == 1
    assert p.sign(x - 3, [x]) == -1
    q = p.extend(CRootOf(x**3 - 2, 0))
    assert q.degree == 6
    assert q.sign(y**3 - 2, [x, y]) == 0 and q.sign(2*x**2 - 9, [x, y]) == 0
    assert q.sign(x*y - 3, [x, y]) == -1
    q = origin.extend(CRootOf(x**3 - 2, 0)).extend(r0)
    assert q.degree == 6
    assert q.sign(2*y**2 - 9, [x, y]) == 0 and q.sign(y, [x, y]) == -1

    # several quadratic extensions
    pp = origin.extend(CRootOf(x**2 - 2, 1)).extend(CRootOf(x**2 - 3, 1))
    pp = pp.extend(CRootOf(x**2 - 5, 1))
    assert pp.degree == 8
    assert pp.sign(x**2*y**2*z**2 - 30, [x, y, z]) == 0
    assert pp.sign(x**2*y**2*z**2 - 29, [x, y, z]) == 1
    assert pp.sign(x*y*z - 6, [x, y, z]) == -1
    assert pp.sign(x*y*z - 5, [x, y, z]) == 1
    vals = [e.evalf(20) for e in pp.as_exprs()]
    assert all(abs(v - sqrt(n).evalf(20)) < 1e-15 for v, n in zip(vals, [2, 3, 5]))
