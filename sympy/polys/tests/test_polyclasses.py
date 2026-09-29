"""Tests for OO layer of several polynomial representations. """
from __future__ import annotations

from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys.domains import ZZ, QQ
from sympy.polys.polyclasses import DMP, SMP, DMF, ANP
from sympy.polys.polytools import Poly
from sympy.polys.polyerrors import (CoercionFailed, ExactQuotientFailed,
                                    NotInvertible, PolynomialError,
                                    UnificationFailed)
from sympy.polys.specialpolys import f_polys
from sympy.testing.pytest import raises, warns_deprecated_sympy

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = [ f.to_dense() for f in f_polys() ]

def test_DMP___init__():
    f = DMP([[ZZ(0)], [], [ZZ(0), ZZ(1), ZZ(2)], [ZZ(3)]], ZZ)

    assert f._rep == [[1, 2], [3]]
    assert f.dom == ZZ
    assert f.lev == 1

    f = DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ, 1)

    assert f._rep == [[1, 2], [3]]
    assert f.dom == ZZ
    assert f.lev == 1

    f = DMP.from_dict({(1, 1): ZZ(1), (0, 0): ZZ(2)}, 1, ZZ)

    assert f._rep == [[1, 0], [2]]
    assert f.dom == ZZ
    assert f.lev == 1


def test_DMP_rep_deprecation():
    f = DMP([1, 2, 3], ZZ)

    with warns_deprecated_sympy():
        assert f.rep == [1, 2, 3]


def test_DMP___eq__():
    assert DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ) == \
        DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ)

    assert DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ) == \
        DMP([[QQ(1), QQ(2)], [QQ(3)]], QQ)
    assert DMP([[QQ(1), QQ(2)], [QQ(3)]], QQ) == \
        DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ)

    assert DMP([[[ZZ(1)]]], ZZ) != DMP([[ZZ(1)]], ZZ)
    assert DMP([[ZZ(1)]], ZZ) != DMP([[[ZZ(1)]]], ZZ)


def test_DMP___bool__():
    assert bool(DMP([[]], ZZ)) is False
    assert bool(DMP([[ZZ(1)]], ZZ)) is True


def test_DMP_to_dict():
    f = DMP([[ZZ(3)], [], [ZZ(2)], [], [ZZ(8)]], ZZ)

    assert f.to_dict() == \
        {(4, 0): 3, (2, 0): 2, (0, 0): 8}
    assert f.to_sympy_dict() == \
        {(4, 0): ZZ.to_sympy(3), (2, 0): ZZ.to_sympy(2), (0, 0):
         ZZ.to_sympy(8)}


def test_DMP_properties():
    assert DMP([[]], ZZ).is_zero is True
    assert DMP([[ZZ(1)]], ZZ).is_zero is False

    assert DMP([[ZZ(1)]], ZZ).is_one is True
    assert DMP([[ZZ(2)]], ZZ).is_one is False

    assert DMP([[ZZ(1)]], ZZ).is_ground is True
    assert DMP([[ZZ(1)], [ZZ(2)], [ZZ(1)]], ZZ).is_ground is False

    assert DMP([[ZZ(1)], [ZZ(2), ZZ(0)], [ZZ(1), ZZ(0)]], ZZ).is_sqf is True
    assert DMP([[ZZ(1)], [ZZ(2), ZZ(0)], [ZZ(1), ZZ(0), ZZ(0)]], ZZ).is_sqf is False

    assert DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ).is_monic is True
    assert DMP([[ZZ(2), ZZ(2)], [ZZ(3)]], ZZ).is_monic is False

    assert DMP([[ZZ(1), ZZ(2)], [ZZ(3)]], ZZ).is_primitive is True
    assert DMP([[ZZ(2), ZZ(4)], [ZZ(6)]], ZZ).is_primitive is False


def test_DMP_arithmetics():
    f = DMP([[ZZ(2)], [ZZ(2), ZZ(0)]], ZZ)

    assert f.mul_ground(2) == DMP([[ZZ(4)], [ZZ(4), ZZ(0)]], ZZ)
    assert f.quo_ground(2) == DMP([[ZZ(1)], [ZZ(1), ZZ(0)]], ZZ)

    raises(ExactQuotientFailed, lambda: f.exquo_ground(3))

    f = DMP([[ZZ(-5)]], ZZ)
    g = DMP([[ZZ(5)]], ZZ)

    assert f.abs() == g
    assert abs(f) == g

    assert g.neg() == f
    assert -g == f

    h = DMP([[]], ZZ)

    assert f.add(g) == h
    assert f + g == h
    assert g + f == h
    assert f + 5 == h
    assert 5 + f == h

    h = DMP([[ZZ(-10)]], ZZ)

    assert f.sub(g) == h
    assert f - g == h
    assert g - f == -h
    assert f - 5 == h
    assert 5 - f == -h

    h = DMP([[ZZ(-25)]], ZZ)

    assert f.mul(g) == h
    assert f * g == h
    assert g * f == h
    assert f * 5 == h
    assert 5 * f == h

    h = DMP([[ZZ(25)]], ZZ)

    assert f.sqr() == h
    assert f.pow(2) == h
    assert f**2 == h

    raises(TypeError, lambda: f.pow('x'))

    f = DMP([[ZZ(1)], [], [ZZ(1), ZZ(0), ZZ(0)]], ZZ)
    g = DMP([[ZZ(2)], [ZZ(-2), ZZ(0)]], ZZ)

    q = DMP([[ZZ(2)], [ZZ(2), ZZ(0)]], ZZ)
    r = DMP([[ZZ(8), ZZ(0), ZZ(0)]], ZZ)

    assert f.pdiv(g) == (q, r)
    assert f.pquo(g) == q
    assert f.prem(g) == r

    raises(ExactQuotientFailed, lambda: f.pexquo(g))

    f = DMP([[ZZ(1)], [], [ZZ(1), ZZ(0), ZZ(0)]], ZZ)
    g = DMP([[ZZ(1)], [ZZ(-1), ZZ(0)]], ZZ)

    q = DMP([[ZZ(1)], [ZZ(1), ZZ(0)]], ZZ)
    r = DMP([[ZZ(2), ZZ(0), ZZ(0)]], ZZ)

    assert f.div(g) == (q, r)
    assert f.quo(g) == q
    assert f.rem(g) == r

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    raises(ExactQuotientFailed, lambda: f.exquo(g))

    f = DMP([ZZ(1), ZZ(0), ZZ(-1)], ZZ)
    g = DMP([ZZ(2), ZZ(-2)], ZZ)

    q = DMP([], ZZ)
    r = f

    pq = DMP([ZZ(2), ZZ(2)], ZZ)
    pr = DMP([], ZZ)

    assert f.div(g) == (q, r)
    assert f.quo(g) == q
    assert f.rem(g) == r

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    raises(ExactQuotientFailed, lambda: f.exquo(g))

    assert f.pdiv(g) == (pq, pr)
    assert f.pquo(g) == pq
    assert f.prem(g) == pr
    assert f.pexquo(g) == pq


def test_DMP_functionality():
    f = DMP([[ZZ(1)], [ZZ(2), ZZ(0)], [ZZ(1), ZZ(0), ZZ(0)]], ZZ)
    g = DMP([[ZZ(1)], [ZZ(1), ZZ(0)]], ZZ)
    h = DMP([[ZZ(1)]], ZZ)

    assert f.degree() == 2
    assert f.degree_list() == (2, 2)
    assert f.total_degree() == 2

    assert f.LC() == ZZ(1)
    assert f.TC() == ZZ(0)
    assert f.nth(1, 1) == ZZ(2)

    raises(TypeError, lambda: f.nth(0, 'x'))

    assert f.max_norm() == 2
    assert f.l1_norm() == 4

    u = DMP([[ZZ(2)], [ZZ(2), ZZ(0)]], ZZ)

    assert f.diff(m=1, j=0) == u
    assert f.diff(m=1, j=1) == u

    raises(TypeError, lambda: f.diff(m='x', j=0))

    u = DMP([ZZ(1), ZZ(2), ZZ(1)], ZZ)
    v = DMP([ZZ(1), ZZ(2), ZZ(1)], ZZ)

    assert f.eval(a=1, j=0) == u
    assert f.eval(a=1, j=1) == v

    assert f.eval(1).eval(1) == ZZ(4)

    assert f.cofactors(g) == (g, g, h)
    assert f.gcd(g) == g
    assert f.lcm(g) == f

    u = DMP([[QQ(45), QQ(30), QQ(5)]], QQ)
    v = DMP([[QQ(1), QQ(2, 3), QQ(1, 9)]], QQ)

    assert u.monic() == v

    assert (4*f).content() == ZZ(4)
    assert (4*f).primitive() == (ZZ(4), f)

    f = DMP([QQ(1,3), QQ(1)], QQ)
    g = DMP([QQ(1,7), QQ(1)], QQ)

    assert f.cancel(g) == f.cancel(g, include=True) == (
        DMP([QQ(7), QQ(21)], QQ),
        DMP([QQ(3), QQ(21)], QQ)
    )
    assert f.cancel(g, include=False) == (
        QQ(7),
        QQ(3),
        DMP([QQ(1), QQ(3)], QQ),
        DMP([QQ(1), QQ(7)], QQ)
    )

    f = DMP([[ZZ(1)], [ZZ(2)], [ZZ(3)], [ZZ(4)], [ZZ(5)], [ZZ(6)]], ZZ)

    assert f.trunc(3) == DMP([[ZZ(1)], [ZZ(-1)], [], [ZZ(1)], [ZZ(-1)], []], ZZ)

    f = DMP(f_4, ZZ)

    assert f.sqf_part() == -f
    assert f.sqf_list() == (ZZ(-1), [(-f, 1)])

    f = DMP([[ZZ(-1)], [], [], [ZZ(5)]], ZZ)
    g = DMP([[ZZ(3), ZZ(1)], [], []], ZZ)
    h = DMP([[ZZ(45), ZZ(30), ZZ(5)]], ZZ)

    r = DMP([ZZ(675), ZZ(675), ZZ(225), ZZ(25)], ZZ)

    assert f.subresultants(g) == [f, g, h]
    assert f.resultant(g) == r

    f = DMP([ZZ(1), ZZ(3), ZZ(9), ZZ(-13)], ZZ)

    assert f.discriminant() == -11664

    f = DMP([QQ(2), QQ(0)], QQ)
    g = DMP([QQ(1), QQ(0), QQ(-16)], QQ)

    s = DMP([QQ(1, 32), QQ(0)], QQ)
    t = DMP([QQ(-1, 16)], QQ)
    h = DMP([QQ(1)], QQ)

    assert f.half_gcdex(g) == (s, h)
    assert f.gcdex(g) == (s, t, h)

    assert f.invert(g) == s

    f = DMP([[QQ(1)], [QQ(2)], [QQ(3)]], QQ)

    raises(ValueError, lambda: f.half_gcdex(f))
    raises(ValueError, lambda: f.gcdex(f))

    raises(ValueError, lambda: f.invert(f))

    f = DMP(ZZ.map([1, 0, 20, 0, 150, 0, 500, 0, 625, -2, 0, -10, 9]), ZZ)
    g = DMP([ZZ(1), ZZ(0), ZZ(0), ZZ(-2), ZZ(9)], ZZ)
    h = DMP([ZZ(1), ZZ(0), ZZ(5), ZZ(0)], ZZ)

    assert g.compose(h) == f
    assert f.decompose() == [g, h]

    f = DMP([[QQ(1)], [QQ(2)], [QQ(3)]], QQ)

    raises(ValueError, lambda: f.decompose())
    raises(ValueError, lambda: f.sturm())


def test_DMP_exclude():
    f = [[[[[[[[[[[[[[[[[[[[[[[[[[ZZ(1)]], [[]]]]]]]]]]]]]]]]]]]]]]]]]]
    J = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 24, 25]

    assert DMP(f, ZZ).exclude() == (J, DMP([ZZ(1), ZZ(0)], ZZ))
    assert DMP([[ZZ(1)], [ZZ(1), ZZ(0)]], ZZ).exclude() ==\
            ([], DMP([[ZZ(1)], [ZZ(1), ZZ(0)]], ZZ))


def test_DMF__init__():
    f = DMF(([[0], [], [0, 1, 2], [3]], [[1, 2, 3]]), ZZ)

    assert f.num == [[1, 2], [3]]
    assert f.den == [[1, 2, 3]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1, 2], [3]], [[1, 2, 3]]), ZZ, 1)

    assert f.num == [[1, 2], [3]]
    assert f.den == [[1, 2, 3]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[-1], [-2]], [[3], [-4]]), ZZ)

    assert f.num == [[-1], [-2]]
    assert f.den == [[3], [-4]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1], [2]], [[-3], [4]]), ZZ)

    assert f.num == [[-1], [-2]]
    assert f.den == [[3], [-4]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1], [2]], [[-3], [4]]), ZZ)

    assert f.num == [[-1], [-2]]
    assert f.den == [[3], [-4]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[]], [[-3], [4]]), ZZ)

    assert f.num == [[]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(17, ZZ, 1)

    assert f.num == [[17]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1], [2]]), ZZ)

    assert f.num == [[1], [2]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF([[0], [], [0, 1, 2], [3]], ZZ)

    assert f.num == [[1, 2], [3]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF({(1, 1): 1, (0, 0): 2}, ZZ, 1)

    assert f.num == [[1, 0], [2]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[QQ(1)], [QQ(2)]], [[-QQ(3)], [QQ(4)]]), QQ)

    assert f.num == [[-QQ(1)], [-QQ(2)]]
    assert f.den == [[QQ(3)], [-QQ(4)]]
    assert f.lev == 1
    assert f.dom == QQ

    f = DMF(([[QQ(1, 5)], [QQ(2, 5)]], [[-QQ(3, 7)], [QQ(4, 7)]]), QQ)

    assert f.num == [[-QQ(7)], [-QQ(14)]]
    assert f.den == [[QQ(15)], [-QQ(20)]]
    assert f.lev == 1
    assert f.dom == QQ

    raises(ValueError, lambda: DMF(([1], [[1]]), ZZ))
    raises(ZeroDivisionError, lambda: DMF(([1], []), ZZ))


def test_DMF__bool__():
    assert bool(DMF([[]], ZZ)) is False
    assert bool(DMF([[1]], ZZ)) is True


def test_DMF_properties():
    assert DMF([[]], ZZ).is_zero is True
    assert DMF([[]], ZZ).is_one is False

    assert DMF([[1]], ZZ).is_zero is False
    assert DMF([[1]], ZZ).is_one is True

    assert DMF(([[1]], [[2]]), ZZ).is_one is False


def test_DMF_arithmetics():
    f = DMF([[7], [-9]], ZZ)
    g = DMF([[-7], [9]], ZZ)

    assert f.neg() == -f == g

    f = DMF(([[1]], [[1], []]), ZZ)
    g = DMF(([[1]], [[1, 0]]), ZZ)

    h = DMF(([[1], [1, 0]], [[1, 0], []]), ZZ)

    assert f.add(g) == f + g == h
    assert g.add(f) == g + f == h

    h = DMF(([[-1], [1, 0]], [[1, 0], []]), ZZ)

    assert f.sub(g) == f - g == h

    h = DMF(([[1]], [[1, 0], []]), ZZ)

    assert f.mul(g) == f*g == h
    assert g.mul(f) == g*f == h

    h = DMF(([[1, 0]], [[1], []]), ZZ)

    assert f.quo(g) == f/g == h

    h = DMF(([[1]], [[1], [], [], []]), ZZ)

    assert f.pow(3) == f**3 == h

    h = DMF(([[1]], [[1, 0, 0, 0]]), ZZ)

    assert g.pow(3) == g**3 == h

    h = DMF(([[1, 0]], [[1]]), ZZ)

    assert g.pow(-1) == g**-1 == h


def test_ANP___init__():
    rep = [QQ(1), QQ(1)]
    mod = [QQ(1), QQ(0), QQ(1)]

    f = ANP(rep, mod, QQ)

    assert f.to_list() == [QQ(1), QQ(1)]
    assert f.mod_to_list() == [QQ(1), QQ(0), QQ(1)]
    assert f.dom == QQ

    rep = {1: QQ(1), 0: QQ(1)}
    mod = {2: QQ(1), 0: QQ(1)}

    f = ANP(rep, mod, QQ)

    assert f.to_list() == [QQ(1), QQ(1)]
    assert f.mod_to_list() == [QQ(1), QQ(0), QQ(1)]
    assert f.dom == QQ

    f = ANP(1, mod, QQ)

    assert f.to_list() == [QQ(1)]
    assert f.mod_to_list() == [QQ(1), QQ(0), QQ(1)]
    assert f.dom == QQ

    f = ANP([1, 0.5], mod, QQ)

    assert all(QQ.of_type(a) for a in f.to_list())

    raises(CoercionFailed, lambda: ANP([sqrt(2)], mod, QQ))

    # https://github.com/sympy/sympy/issues/28182
    f = ANP({(0,): -2, (2,): 1}, [1, 0, -2], QQ)
    assert f.to_list() == []


def test_ANP___eq__():
    a = ANP([QQ(1), QQ(1)], [QQ(1), QQ(0), QQ(1)], QQ)
    b = ANP([QQ(1), QQ(1)], [QQ(1), QQ(0), QQ(2)], QQ)

    assert (a == a) is True
    assert (a != a) is False

    assert (a == b) is False
    assert (a != b) is True

    b = ANP([QQ(1), QQ(2)], [QQ(1), QQ(0), QQ(1)], QQ)

    assert (a == b) is False
    assert (a != b) is True


def test_ANP___bool__():
    assert bool(ANP([], [QQ(1), QQ(0), QQ(1)], QQ)) is False
    assert bool(ANP([QQ(1)], [QQ(1), QQ(0), QQ(1)], QQ)) is True


def test_ANP_properties():
    mod = [QQ(1), QQ(0), QQ(1)]

    assert ANP([QQ(0)], mod, QQ).is_zero is True
    assert ANP([QQ(1)], mod, QQ).is_zero is False

    assert ANP([QQ(1)], mod, QQ).is_one is True
    assert ANP([QQ(2)], mod, QQ).is_one is False


def test_ANP_arithmetics():
    mod = [QQ(1), QQ(0), QQ(0), QQ(-2)]

    a = ANP([QQ(2), QQ(-1), QQ(1)], mod, QQ)
    b = ANP([QQ(1), QQ(2)], mod, QQ)

    c = ANP([QQ(-2), QQ(1), QQ(-1)], mod, QQ)

    assert a.neg() == -a == c

    c = ANP([QQ(2), QQ(0), QQ(3)], mod, QQ)

    assert a.add(b) == a + b == c
    assert b.add(a) == b + a == c

    c = ANP([QQ(2), QQ(-2), QQ(-1)], mod, QQ)

    assert a.sub(b) == a - b == c

    c = ANP([QQ(-2), QQ(2), QQ(1)], mod, QQ)

    assert b.sub(a) == b - a == c

    c = ANP([QQ(3), QQ(-1), QQ(6)], mod, QQ)

    assert a.mul(b) == a*b == c
    assert b.mul(a) == b*a == c

    c = ANP([QQ(-1, 43), QQ(9, 43), QQ(5, 43)], mod, QQ)

    assert a.pow(0) == a**(0) == ANP(1, mod, QQ)
    assert a.pow(1) == a**(1) == a

    assert a.pow(-1) == a**(-1) == c

    assert a.quo(a) == a.mul(a.pow(-1)) == a*a**(-1) == ANP(1, mod, QQ)

    c = ANP([], [1, 0, 0, -2], QQ)
    r1 = a.rem(b)

    (q, r2) = a.div(b)

    assert r1 == r2 == c == a % b

    raises(NotInvertible, lambda: a.div(c))
    raises(NotInvertible, lambda: a.rem(c))

    # Comparison with "hard-coded" value fails despite looking identical
    # from sympy import Rational
    # c = ANP([Rational(11, 10), Rational(-1, 5), Rational(-3, 5)], [1, 0, 0, -2], QQ)

    assert q == a/b # == c

def test_ANP_unify():
    mod_z = [ZZ(1), ZZ(0), ZZ(-2)]
    mod_q = [QQ(1), QQ(0), QQ(-2)]

    a = ANP([QQ(1)], mod_q, QQ)
    b = ANP([ZZ(1)], mod_z, ZZ)

    assert a.unify(b)[0] == QQ
    assert b.unify(a)[0] == QQ
    assert a.unify(a)[0] == QQ
    assert b.unify(b)[0] == ZZ

    assert a.unify_ANP(b)[-1] == QQ
    assert b.unify_ANP(a)[-1] == QQ
    assert a.unify_ANP(a)[-1] == QQ
    assert b.unify_ANP(b)[-1] == ZZ


def test_zero_poly():
    from sympy import Symbol
    x = Symbol('x')

    R_old = ZZ.old_poly_ring(x)
    zero_poly_old = R_old(0)
    cont_old, prim_old = zero_poly_old.primitive()

    assert cont_old == 0
    assert prim_old == zero_poly_old
    assert prim_old.is_primitive is False




def test_SMP_constructors():
    f = SMP({(2, 0): 1, (0, 1): 2}, ZZ)

    assert f.dom == ZZ
    assert f.lev == 1
    assert f.to_dict() == {(2, 0): 1, (0, 1): 2}

    f = SMP.from_dict({(2, 0): 1, (0, 1): 2}, 1, ZZ)

    assert f.dom == ZZ
    assert f.lev == 1
    assert f.to_dict() == {(2, 0): 1, (0, 1): 2}

    assert SMP.zero(1, ZZ).to_dict() == {}
    assert SMP.one(1, ZZ).to_dict() == {(0, 0): 1}

    raises(ValueError, lambda: SMP({}, ZZ))


def test_SMP_from_poly():
    x, y = symbols('x y')

    p = Poly(x**2 + y, x, y)
    f = SMP.from_poly(p)

    assert f.dom == p.domain
    assert f.lev == p.rep.lev
    assert f.to_dict() == p.rep.to_dict()

    p = Poly.new(f, x, y)
    assert SMP.from_poly(p) is f


def test_SMP_from_expr():
    x, y, z = symbols('x y z')

    expr = x**2 + 2*x - 5*y - 7*y**3*x**3 + 3
    f, gens = SMP.from_expr(expr)

    assert gens == (x, y)
    assert f.dom == ZZ
    assert f.to_dict() == {
        (3, 3): -7,
        (2, 0): 1,
        (1, 0): 2,
        (0, 1): -5,
        (0, 0): 3,
    }

    f, gens = SMP.from_expr(x + y, x)

    assert gens == (x,)
    assert f.to_sympy_dict() == {(1,): 1, (0,): y}
    assert str(f.dom) == 'ZZ[y]'


def test_SMP_constructor_errors_and_domain():
    x = symbols('x')

    raises(PolynomialError, lambda: SMP.from_dict({(1, 0): 1}, 0, ZZ))
    raises(PolynomialError, lambda: SMP.from_expr(1))

    f, gens = SMP.from_expr(x/2 + 1, x, domain=QQ)
    assert gens == (x,)
    assert f.dom == QQ
    assert f.to_dict() == {(1,): QQ(1, 2), (0,): QQ.one}


def test_SMP_equality_and_unification():
    f = SMP.from_dict({(1,): 1, (0,): 2}, 0, ZZ)
    g = SMP.from_dict({(1,): QQ.one, (0,): QQ(2)}, 0, QQ)

    assert f == g
    assert f.eq(g)
    assert not f.eq(g, strict=True)

    F, G = f.unify_SMP(g)
    assert F.dom == G.dom == QQ
    assert F.to_dict() == G.to_dict()

    h = SMP.from_dict({(1, 0): 1}, 1, ZZ)
    assert f != h
    raises(UnificationFailed, lambda: f.unify_SMP(h))


def test_SMP_arithmetic_and_terms():
    f = SMP.from_dict({(2,): 1, (0,): 1}, 0, ZZ)
    g = SMP.from_dict({(1,): 1, (0,): -1}, 0, ZZ)

    assert f.add_ground(2).to_dict() == {(2,): 1, (0,): 3}
    assert f.sub_ground(1).to_dict() == {(2,): 1}
    assert f.mul_ground(3).to_dict() == {(2,): 3, (0,): 3}
    assert f.neg().to_dict() == {(2,): -1, (0,): -1}

    assert f.add(g).to_dict() == {(2,): 1, (1,): 1}
    assert f.sub(g).to_dict() == {(2,): 1, (1,): -1, (0,): 2}
    assert f.mul(g).to_dict() == {(3,): 1, (2,): -1, (1,): 1, (0,): -1}
    assert g.sqr().to_dict() == {(2,): 1, (1,): -2, (0,): 1}
    assert g.pow(3).to_dict() == {(3,): 1, (2,): -3, (1,): 3, (0,): -1}

    raises(TypeError, lambda: g.pow(QQ(2)))
    raises(ValueError, lambda: g.pow(-1))

    assert f.coeffs() == [ZZ.one, ZZ.one]
    assert f.monoms() == [(2,), (0,)]
    assert f.terms() == [((2,), ZZ.one), ((0,), ZZ.one)]

    z = SMP.zero(0, ZZ)
    assert z.terms() == [((0,), ZZ.zero)]
    assert z.all_coeffs() == [ZZ.zero]
    assert z.all_monoms() == [(0,)]
    assert z.all_terms() == [((0,), ZZ.zero)]

    h = SMP.from_dict({(2, 1): 1}, 1, ZZ)
    raises(PolynomialError, lambda: h.all_coeffs())
    raises(PolynomialError, lambda: h.all_monoms())


def test_SMP_degrees_coefficients_and_conversion():
    f = SMP.from_dict({(3, 1): 2, (0, 4): 3, (0, 0): 6}, 1, ZZ)

    assert f.degree(0) == 3
    assert f.degree(1) == 4
    assert f.degree_list() == (3, 4)
    assert f.total_degree() == 4
    assert f.LC() == ZZ(2)
    assert f.TC() == ZZ(6)
    assert f.nth(3, 1) == ZZ(2)
    assert f.nth(1, 1) == ZZ.zero

    raises(TypeError, lambda: f.degree(QQ.one))
    raises(TypeError, lambda: f.nth(1, QQ.one))

    assert f.to_list() == DMP.from_dict(f.to_dict(), f.lev, f.dom).to_list()
    assert f.to_tuple() == DMP.from_dict(f.to_dict(), f.lev, f.dom).to_tuple()
    assert f.to_dict(zero=True) == f.to_dict()

    z = SMP.zero(1, ZZ)
    assert z.to_dict(zero=True) == {(0, 0): ZZ.zero}

    q = SMP.from_dict({(1,): QQ(1, 2), (0,): QQ(1, 3)}, 0, QQ)
    c, q0 = q.clear_denoms()
    assert c == ZZ(6)
    assert q0.to_dict() == {(1,): QQ(3), (0,): QQ(2)}

    cont, prim = SMP.from_dict({(1,): 6, (0,): 9}, 0, ZZ).primitive()
    assert cont == ZZ(3)
    assert prim.to_dict() == {(1,): 2, (0,): 3}


def test_SMP_calculus_eval_and_trunc():
    f = SMP.from_dict({(3,): QQ(2), (1,): QQ(5), (0,): QQ(7)}, 0, QQ)

    assert f.diff().to_dict() == {(2,): QQ(6), (0,): QQ(5)}
    assert f.integrate().to_dict() == {
        (4,): QQ(1, 2), (2,): QQ(5, 2), (1,): QQ(7)}

    assert f.eval(QQ(2)) == QQ(33)

    g = SMP.from_dict({(1, 1): 2, (0, 1): 3, (0, 0): 1}, 1, ZZ)
    assert g.eval(2, 0).to_dict() == {(1,): 7, (0,): 1}

    h = SMP.from_dict({(2,): 7, (1,): -6, (0,): 10}, 0, ZZ)
    assert h.trunc(5).to_dict() == {(2,): 2, (1,): -1}


def test_SMP_gcd_lcm_and_polynomial_algorithms():
    x, y = symbols('x y')

    f = Poly(x**3 - x, x)
    g = Poly(x**2 - 1, x)
    F = SMP.from_poly(f)
    G = SMP.from_poly(g)

    assert F.gcd(G).to_dict() == f.rep.gcd(g.rep).to_dict()
    assert F.lcm(G).to_dict() == f.rep.lcm(g.rep).to_dict()
    assert F.lcm(SMP.zero(0, ZZ)).is_zero

    Fq = F.convert(QQ)
    Gq = G.convert(QQ)
    assert Fq.lcm(Gq).to_dict() == f.rep.convert(QQ).lcm(
        g.rep.convert(QQ)).to_dict()

    assert [h.to_dict() for h in F.subresultants(G)] == [
        h.to_dict() for h in f.rep.subresultants(g.rep)]

    r = F.resultant(G)
    assert r == f.rep.resultant(g.rep)

    r, prs = F.resultant(G, includePRS=True)
    dr, dprs = f.rep.resultant(g.rep, includePRS=True)
    assert r == dr
    assert [h.to_dict() for h in prs] == [h.to_dict() for h in dprs]

    fm = Poly(x**2 + y*x + 1, x, y)
    gm = Poly(x + y, x, y)
    Fm = SMP.from_poly(fm)
    Gm = SMP.from_poly(gm)
    assert Fm.resultant(Gm).to_dict() == fm.rep.resultant(gm.rep).to_dict()

    d = SMP.from_poly(Poly(x**3 - 2*x + 1, x))
    assert d.discriminant() == Poly(x**3 - 2*x + 1, x).rep.discriminant()

    sf = SMP.from_poly(Poly((x - 1)**2*(x + 2)**3, x))
    coeff, factors = sf.sqf_list()
    dcoeff, dfactors = Poly((x - 1)**2*(x + 2)**3, x).rep.sqf_list()
    assert coeff == dcoeff
    assert [(h.to_dict(), k) for h, k in factors] == [
        (h.to_dict(), k) for h, k in dfactors]

    sfm = SMP.from_poly(Poly((x + y)**2*(x - y), x, y))
    coeff, factors = sfm.sqf_list()
    dcoeff, dfactors = Poly((x + y)**2*(x - y), x, y).rep.sqf_list()
    assert coeff == dcoeff
    assert [(h.to_dict(), k) for h, k in factors] == [
        (h.to_dict(), k) for h, k in dfactors]

    coeff, factors = sf.factor_list()
    dcoeff, dfactors = Poly((x - 1)**2*(x + 2)**3, x).rep.factor_list()
    assert coeff == dcoeff
    assert [(h.to_dict(), k) for h, k in factors] == [
        (h.to_dict(), k) for h, k in dfactors]


def test_SMP_properties():
    z = SMP.zero(0, ZZ)
    one = SMP.one(0, ZZ)
    f = SMP.from_dict({(1,): 1, (0,): 2}, 0, ZZ)
    g = SMP.from_dict({(2,): 1, (0,): 1}, 0, ZZ)
    h = SMP.from_dict({(1,): 2, (0,): 4}, 0, ZZ)

    assert z.is_zero
    assert not f.is_zero
    assert one.is_one
    assert not f.is_one
    assert one.is_ground
    assert not f.is_ground
    assert f.is_monic
    assert f.is_primitive
    assert not h.is_primitive
    assert f.is_linear
    assert not g.is_linear
    assert g.is_quadratic
    assert f.is_monomial is False
    assert SMP.from_dict({(3,): 2}, 0, ZZ).is_monomial
