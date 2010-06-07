"""Tests for OO layer of several polynomial representations. """

from sympy.polys.polyclasses import (
    GFP, init_normal_GFP,
    DMP, init_normal_DMP,
    SDP, init_normal_SDP,
    DMF, init_normal_DMF,
    ANP, init_normal_ANP,
)

from sympy.polys.algebratools import ZZ, QQ
from sympy.polys.specialpolys import f_4

from sympy.polys.polyerrors import (
    ExactQuotientFailed,
)

from sympy import raises

def test_DMP___init__():
    f = DMP([[0],[],[0,1,2],[3]], ZZ)

    assert f.rep == [[1,2],[3]]
    assert f.dom == ZZ
    assert f.lev == 1

    f = DMP([[1,2],[3]], ZZ, 1)

    assert f.rep == [[1,2],[3]]
    assert f.dom == ZZ
    assert f.lev == 1

    f = DMP({(1,1): 1, (0,0): 2}, ZZ, 1)

    assert f.rep == [[1,0],[2]]
    assert f.dom == ZZ
    assert f.lev == 1

def test_DMP___eq__():
    assert DMP([[ZZ(1),ZZ(2)],[ZZ(3)]], ZZ) == \
           DMP([[ZZ(1),ZZ(2)],[ZZ(3)]], ZZ)

    assert DMP([[ZZ(1),ZZ(2)],[ZZ(3)]], ZZ) == \
           DMP([[QQ(1),QQ(2)],[QQ(3)]], QQ)
    assert DMP([[QQ(1),QQ(2)],[QQ(3)]], QQ) == \
           DMP([[ZZ(1),ZZ(2)],[ZZ(3)]], ZZ)

    assert DMP([[[ZZ(1)]]], ZZ) != DMP([[ZZ(1)]], ZZ)
    assert DMP([[ZZ(1)]], ZZ) != DMP([[[ZZ(1)]]], ZZ)

def test_DMP___bool__():
    assert bool(DMP([[]], ZZ)) == False
    assert bool(DMP([[1]], ZZ)) == True

def test_DMP_to_dict():
    f = DMP([[3],[],[2],[],[8]], ZZ)

    assert f.to_dict() == \
        {(4, 0): 3, (2, 0): 2, (0, 0): 8}
    assert f.to_sympy_dict() == \
        {(4, 0): ZZ.to_sympy(3), (2, 0): ZZ.to_sympy(2), (0, 0): ZZ.to_sympy(8)}

def test_DMP_properties():
    assert DMP([[]], ZZ).is_zero == True
    assert DMP([[1]], ZZ).is_zero == False

    assert DMP([[1]], ZZ).is_one == True
    assert DMP([[2]], ZZ).is_one == False

    assert DMP([[1]], ZZ).is_ground == True
    assert DMP([[1],[2],[1]], ZZ).is_ground == False

    assert DMP([[1],[2,0],[1,0]], ZZ).is_sqf == True
    assert DMP([[1],[2,0],[1,0,0]], ZZ).is_sqf == False

    assert DMP([[1,2],[3]], ZZ).is_monic == True
    assert DMP([[2,2],[3]], ZZ).is_monic == False

    assert DMP([[1,2],[3]], ZZ).is_primitive == True
    assert DMP([[2,4],[6]], ZZ).is_primitive == False

def test_DMP_arithmetics():
    f = DMP([[2],[2,0]], ZZ)

    assert f.mul_ground(2) == DMP([[4],[4,0]], ZZ)
    assert f.exquo_ground(2) == DMP([[1],[1,0]], ZZ)

    raises(ExactQuotientFailed, 'f.quo_ground(3)')

    f = DMP([[-5]], ZZ)
    g = DMP([[5]], ZZ)

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

    h = DMP([[-10]], ZZ)

    assert f.sub(g) == h
    assert f - g ==  h
    assert g - f == -h
    assert f - 5 ==  h
    assert 5 - f == -h

    h = DMP([[-25]], ZZ)

    assert f.mul(g) == h
    assert f * g == h
    assert g * f == h
    assert f * 5 == h
    assert 5 * f == h

    h = DMP([[25]], ZZ)

    assert f.sqr() == h
    assert f.pow(2) == h
    assert f**2 == h

    raises(TypeError, "f.pow('x')")

    f = DMP([[1],[],[1,0,0]], ZZ)
    g = DMP([[2],[-2,0]], ZZ)

    q = DMP([[2],[2,0]], ZZ)
    r = DMP([[8,0,0]], ZZ)

    assert f.pdiv(g) == (q, r)
    assert f.pexquo(g) == q
    assert f.prem(g) == r

    raises(ExactQuotientFailed, 'f.pquo(g)')

    f = DMP([[1],[],[1,0,0]], ZZ)
    g = DMP([[1],[-1,0]], ZZ)

    q = DMP([[1],[1,0]], ZZ)
    r = DMP([[2,0,0]], ZZ)

    assert f.div(g) == (q, r)
    assert f.exquo(g) == q
    assert f.rem(g) == r

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    raises(ExactQuotientFailed, 'f.quo(g)')

def test_DMP_functionality():
    f = DMP([[1],[2,0],[1,0,0]], ZZ)
    g = DMP([[1],[1,0]], ZZ)
    h = DMP([[1]], ZZ)

    assert f.degree() == 2
    assert f.degree_list() == (2, 2)
    assert f.total_degree() == 4

    assert f.LC() == ZZ(1)
    assert f.TC() == ZZ(0)
    assert f.nth(1, 1) == ZZ(2)

    raises(TypeError, "f.nth(0, 'x')")

    assert f.max_norm() == 2
    assert f.l1_norm() == 4

    u = DMP([[2],[2,0]], ZZ)

    assert f.diff(m=1, j=0) == u
    assert f.diff(m=1, j=1) == u

    raises(TypeError, "f.diff(m='x', j=0)")

    u = DMP([1,2,1], ZZ)
    v = DMP([1,2,1], ZZ)

    assert f.eval(a=1, j=0) == u
    assert f.eval(a=1, j=1) == v

    assert f.eval(1).eval(1) == ZZ(4)

    assert f.cofactors(g) == (g, g, h)
    assert f.gcd(g) == g
    assert f.lcm(g) == f

    u = DMP([[QQ(45),QQ(30),QQ(5)]], QQ)
    v = DMP([[QQ(1),QQ(2,3),QQ(1,9)]], QQ)

    assert u.monic() == v

    assert (4*f).content() == ZZ(4)
    assert (4*f).primitive() == (ZZ(4), f)

    f = DMP([[1],[2],[3],[4],[5],[6]], ZZ)

    assert f.trunc(3) == DMP([[1],[-1],[],[1],[-1],[]], ZZ)

    f = DMP(f_4, ZZ)

    assert f.sqf_part() == -f
    assert f.sqf_list() == (ZZ(-1), [(-f, 1)])

    f = DMP([[-1],[],[],[5]], ZZ)
    g = DMP([[3,1],[],[]], ZZ)
    h = DMP([[45,30,5]], ZZ)

    r = DMP([675,675,225,25], ZZ)

    assert f.subresultants(g) == [f, g, h]
    assert f.resultant(g) == r

    f = DMP([1,3,9,-13], ZZ)

    assert f.discriminant() == -11664

    f = DMP([QQ(2),QQ(0)], QQ)
    g = DMP([QQ(1),QQ(0),QQ(-16)], QQ)

    s = DMP([QQ(1,32),QQ(0)], QQ)
    t = DMP([QQ(-1,16)], QQ)
    h = DMP([QQ(1)], QQ)

    assert f.half_gcdex(g) == (s, h)
    assert f.gcdex(g) == (s, t, h)

    assert f.invert(g) == s

    f = DMP([[1],[2],[3]], QQ)

    raises(ValueError, "f.half_gcdex(f)")
    raises(ValueError, "f.gcdex(f)")

    raises(ValueError, "f.invert(f)")

    f = DMP([1,0,20,0,150,0,500,0,625,-2,0,-10,9], ZZ)
    g = DMP([1,0,0,-2,9], ZZ)
    h = DMP([1,0,5,0], ZZ)

    assert g.compose(h) == f
    assert f.decompose() == [g, h]

    f = DMP([[1],[2],[3]], QQ)

    raises(ValueError, "f.decompose()")
    raises(ValueError, "f.sturm()")

def test_DMP_exclude():
    assert DMP([[[[[[[[[[[[[[[[[[[[[[[[[[1]], [[]]]]]]]]]]]]]]]]]]]]]]]]]],
    ZZ).exclude() == \
        ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 24, 25], DMP([1, 0], ZZ))

    assert DMP([[1], [1, 0]], ZZ).exclude() == ([], DMP([[1], [1, 0]], ZZ))


def test_DMF__init__():
    f = DMF(([[0],[],[0,1,2],[3]], [[1,2,3]]), ZZ)

    assert f.num == [[1,2],[3]]
    assert f.den == [[1,2,3]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1,2],[3]], [[1,2,3]]), ZZ, 1)

    assert f.num == [[1,2],[3]]
    assert f.den == [[1,2,3]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[-1],[-2]],[[3],[-4]]), ZZ)

    assert f.num == [[-1],[-2]]
    assert f.den == [[3],[-4]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1],[2]],[[-3],[4]]), ZZ)

    assert f.num == [[-1],[-2]]
    assert f.den == [[3],[-4]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1],[2]],[[-3],[4]]), ZZ)

    assert f.num == [[-1],[-2]]
    assert f.den == [[3],[-4]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[]],[[-3],[4]]), ZZ)

    assert f.num == [[]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(17, ZZ, 1)

    assert f.num == [[17]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF(([[1],[2]]), ZZ)

    assert f.num == [[1],[2]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF([[0],[],[0,1,2],[3]], ZZ)

    assert f.num == [[1,2],[3]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    f = DMF({(1,1): 1, (0,0): 2}, ZZ, 1)

    assert f.num == [[1,0],[2]]
    assert f.den == [[1]]
    assert f.lev == 1
    assert f.dom == ZZ

    raises(ValueError, "DMF(([1], [[1]]), ZZ)")
    raises(ZeroDivisionError, "DMF(([1], []), ZZ)")

def test_DMF__eq__():
    pass

def test_DMF__bool__():
    assert bool(DMF([[]], ZZ)) == False
    assert bool(DMF([[1]], ZZ)) == True

def test_DMF_properties():
    assert DMF([[]], ZZ).is_zero == True
    assert DMF([[]], ZZ).is_one == False

    assert DMF([[1]], ZZ).is_zero == False
    assert DMF([[1]], ZZ).is_one == True

    assert DMF(([[1]], [[2]]), ZZ).is_one == False

def test_DMF_arithmetics():
    f = DMF([[7],[-9]], ZZ)
    g = DMF([[-7],[9]], ZZ)

    assert f.neg() == -f == g

    f = DMF(([[1]], [[1],[]]), ZZ)
    g = DMF(([[1]], [[1,0]]), ZZ)

    h = DMF(([[1],[1,0]], [[1,0],[]]), ZZ)

    assert f.add(g) == f + g == h
    assert g.add(f) == g + f == h

    h = DMF(([[-1],[1,0]], [[1,0],[]]), ZZ)

    assert f.sub(g) == f - g == h

    h = DMF(([[1]], [[1,0],[]]), ZZ)

    assert f.mul(g) == f*g == h
    assert g.mul(f) == g*f == h

    h = DMF(([[1,0]], [[1],[]]), ZZ)

    assert f.quo(g) == f/g == h

    h = DMF(([[1]], [[1],[],[],[]]), ZZ)

    assert f.pow(3) == f**3 == h

    h = DMF(([[1]], [[1,0,0,0]]), ZZ)

    assert g.pow(3) == g**3 == h

def test_ANP___init__():
    rep = [QQ(1),QQ(1)]
    mod = [QQ(1),QQ(0),QQ(1)]

    f = ANP(rep, mod, QQ)

    assert f.rep == [QQ(1),QQ(1)]
    assert f.mod == [QQ(1),QQ(0),QQ(1)]
    assert f.dom == QQ

    rep = {1: QQ(1), 0: QQ(1)}
    mod = {2: QQ(1), 0: QQ(1)}

    f = ANP(rep, mod, QQ)

    assert f.rep == [QQ(1),QQ(1)]
    assert f.mod == [QQ(1),QQ(0),QQ(1)]
    assert f.dom == QQ

    f = ANP(1, mod, QQ)

    assert f.rep == [QQ(1)]
    assert f.mod == [QQ(1),QQ(0),QQ(1)]
    assert f.dom == QQ

def test_ANP___eq__():
    a = ANP([QQ(1), QQ(1)], [QQ(1),QQ(0),QQ(1)], QQ)
    b = ANP([QQ(1), QQ(1)], [QQ(1),QQ(0),QQ(2)], QQ)

    assert (a == a) == True
    assert (a != a) == False

    assert (a == b) == False
    assert (a != b) == True

    b = ANP([QQ(1), QQ(2)], [QQ(1),QQ(0),QQ(1)], QQ)

    assert (a == b) == False
    assert (a != b) == True

def test_ANP___bool__():
    assert bool(ANP([], [QQ(1),QQ(0),QQ(1)], QQ)) == False
    assert bool(ANP([QQ(1)], [QQ(1),QQ(0),QQ(1)], QQ)) == True

def test_ANP_properties():
    mod = [QQ(1),QQ(0),QQ(1)]

    assert ANP([QQ(0)], mod, QQ).is_zero == True
    assert ANP([QQ(1)], mod, QQ).is_zero == False

    assert ANP([QQ(1)], mod, QQ).is_one == True
    assert ANP([QQ(2)], mod, QQ).is_one == False

def test_ANP_arithmetics():
    mod = [QQ(1),QQ(0),QQ(0),QQ(-2)]

    a = ANP([QQ(2),QQ(-1),QQ(1)], mod, QQ)
    b = ANP([QQ(1),QQ(2)], mod, QQ)

    c = ANP([QQ(-2), QQ(1), QQ(-1)], mod, QQ)

    assert a.neg() == -a == c

    c = ANP([QQ(2), QQ(0), QQ(3)], mod, QQ)

    assert a.add(b) == a+b == c
    assert b.add(a) == b+a == c

    c = ANP([QQ(2), QQ(-2), QQ(-1)], mod, QQ)

    assert a.sub(b) == a-b == c

    c = ANP([QQ(-2), QQ(2), QQ(1)], mod, QQ)

    assert b.sub(a) == b-a == c

    c = ANP([QQ(3), QQ(-1), QQ(6)], mod, QQ)

    assert a.mul(b) == a*b == c
    assert b.mul(a) == b*a == c

    c = ANP([QQ(-1,43), QQ(9,43), QQ(5,43)], mod, QQ)

    assert a.pow(0) == a**(0) == ANP(1, mod, QQ)
    assert a.pow(1) == a**(1) == a

    assert a.pow(-1) == a**(-1) == c

    assert a.quo(a) == a.mul(a.pow(-1)) == a*a**(-1) == ANP(1, mod, QQ)

