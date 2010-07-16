"""Tests for tools for real and complex root isolation and refinement. """

from sympy.polys.rootisolation import (
    dup_refine_real_root,
    dup_isolate_real_roots,
    dup_isolate_real_roots_sqf,
    dup_isolate_real_roots_list,
    dup_count_real_roots,
    dup_count_complex_roots,
    dup_isolate_complex_roots_sqf,
    dup_isolate_all_roots,
    dup_isolate_all_roots_sqf,
)

from sympy.polys.densebasic import (
    dup_from_raw_dict,
)

from sympy.polys.densetools import (
    dup_sqf_part,
)

from sympy.polys.polyerrors import (
    RefinementFailed,
    DomainError,
)

from sympy.polys.domains import ZZ, QQ, EX

from sympy import raises

def test_dup_refine_real_root():
    f = [1,0,-2]

    assert dup_refine_real_root(f, QQ(1), QQ(1), ZZ, steps=1) == (QQ(1), QQ(1))
    assert dup_refine_real_root(f, QQ(1), QQ(1), ZZ, steps=9) == (QQ(1), QQ(1))

    raises(ValueError, "dup_refine_real_root(f, QQ(-2), QQ(2), ZZ)")

    s, t = QQ(1,1), QQ(2,1)

    assert dup_refine_real_root(f, s, t, ZZ, steps=0) == (QQ(1, 1), QQ(2, 1))
    assert dup_refine_real_root(f, s, t, ZZ, steps=1) == (QQ(1, 1), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=2) == (QQ(4, 3), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=3) == (QQ(7, 5), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=4) == (QQ(7, 5), QQ(10, 7))

    s, t = QQ(1,1), QQ(3,2)

    assert dup_refine_real_root(f, s, t, ZZ, steps=0) == (QQ(1, 1), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=1) == (QQ(4, 3), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=2) == (QQ(7, 5), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=3) == (QQ(7, 5), QQ(10, 7))
    assert dup_refine_real_root(f, s, t, ZZ, steps=4) == (QQ(7, 5), QQ(17, 12))

    s, t = QQ(1,1), QQ(5,3)

    assert dup_refine_real_root(f, s, t, ZZ, steps=0) == (QQ(1, 1), QQ(5, 3))
    assert dup_refine_real_root(f, s, t, ZZ, steps=1) == (QQ(1, 1), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=2) == (QQ(7, 5), QQ(3, 2))
    assert dup_refine_real_root(f, s, t, ZZ, steps=3) == (QQ(7, 5), QQ(13, 9))
    assert dup_refine_real_root(f, s, t, ZZ, steps=4) == (QQ(7, 5), QQ(10, 7))

    s, t = QQ(-1,1), QQ(-2,1)

    assert dup_refine_real_root(f, s, t, ZZ, steps=0) == (-QQ(2, 1), -QQ(1, 1))
    assert dup_refine_real_root(f, s, t, ZZ, steps=1) == (-QQ(3, 2), -QQ(1, 1))
    assert dup_refine_real_root(f, s, t, ZZ, steps=2) == (-QQ(3, 2), -QQ(4, 3))
    assert dup_refine_real_root(f, s, t, ZZ, steps=3) == (-QQ(3, 2), -QQ(7, 5))
    assert dup_refine_real_root(f, s, t, ZZ, steps=4) == (-QQ(10, 7), -QQ(7, 5))

    raises(RefinementFailed, "dup_refine_real_root(f, QQ(0), QQ(1), ZZ)")

    s, t, u, v, w = QQ(1), QQ(2), QQ(24,17), QQ(17,12), QQ(7,5)

    assert dup_refine_real_root(f, s, t, ZZ, eps=QQ(1,100)) == (u, v)
    assert dup_refine_real_root(f, s, t, ZZ, steps=6) == (u, v)

    assert dup_refine_real_root(f, s, t, ZZ, eps=QQ(1,100), steps=5) == (w, v)
    assert dup_refine_real_root(f, s, t, ZZ, eps=QQ(1,100), steps=6) == (u, v)
    assert dup_refine_real_root(f, s, t, ZZ, eps=QQ(1,100), steps=7) == (u, v)

    s, t, u, v = QQ(-2), QQ(-1), QQ(-3,2), QQ(-4,3)

    assert dup_refine_real_root([1,0,-2], s, t, ZZ, disjoint=QQ(-5)) == (s, t)
    assert dup_refine_real_root([1,0,-2], s, t, ZZ, disjoint=-v) == (s, t)
    assert dup_refine_real_root([1,0,-2], s, t, ZZ, disjoint=v) == (u, v)

    s, t, u, v = QQ(1), QQ(2), QQ(4,3), QQ(3,2)

    assert dup_refine_real_root([1,0,-2], s, t, ZZ, disjoint=QQ(5)) == (s, t)
    assert dup_refine_real_root([1,0,-2], s, t, ZZ, disjoint=-u) == (s, t)
    assert dup_refine_real_root([1,0,-2], s, t, ZZ, disjoint=u) == (u, v)

def test_dup_isolate_real_roots_sqf():
    assert dup_isolate_real_roots_sqf([], ZZ) == []
    assert dup_isolate_real_roots_sqf([5], ZZ) == []

    assert dup_isolate_real_roots_sqf([1, 1,0], ZZ) == [(-QQ(1), -QQ(1)), (QQ(0), QQ(0))]
    assert dup_isolate_real_roots_sqf([1,-1,0], ZZ) == [( QQ(0),  QQ(0)), (QQ(1), QQ(1))]

    assert dup_isolate_real_roots_sqf([1,0,0,1,1], ZZ) == []

    I = [ (-QQ(2), -QQ(1)), (QQ(1), QQ(2))]

    assert dup_isolate_real_roots_sqf([1,0,-2], ZZ) == I
    assert dup_isolate_real_roots_sqf([-1,0,2], ZZ) == I

    assert dup_isolate_real_roots_sqf([1,-1], ZZ) == \
        [(QQ(1), QQ(1))]
    assert dup_isolate_real_roots_sqf([1,-3,2], ZZ) == \
        [(QQ(1), QQ(1)), (QQ(2), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,-6,11,-6], ZZ) == \
        [(QQ(1), QQ(1)), (QQ(2), QQ(2)), (QQ(3), QQ(3))]
    assert dup_isolate_real_roots_sqf([1,-10,35,-50,24], ZZ) == \
        [(QQ(1), QQ(1)), (QQ(2), QQ(2)), (QQ(3), QQ(3)), (QQ(4), QQ(4))]
    assert dup_isolate_real_roots_sqf([1,-15,85,-225,274,-120], ZZ) == \
        [(QQ(1), QQ(1)), (QQ(2), QQ(2)), (QQ(3), QQ(3)), (QQ(4), QQ(4)), (QQ(5), QQ(5))]

    assert dup_isolate_real_roots_sqf([1,-10], ZZ) == \
        [(QQ(10), QQ(10))]
    assert dup_isolate_real_roots_sqf([1,-30,200], ZZ) == \
        [(QQ(10), QQ(10)), (QQ(20), QQ(20))]
    assert dup_isolate_real_roots_sqf([1,-60,1100,-6000], ZZ) == \
        [(QQ(10), QQ(10)), (QQ(20), QQ(20)), (QQ(30), QQ(30))]
    assert dup_isolate_real_roots_sqf([1,-100,3500,-50000,240000], ZZ) == \
        [(QQ(10), QQ(10)), (QQ(20), QQ(20)), (QQ(30), QQ(30)), (QQ(40), QQ(40))]
    assert dup_isolate_real_roots_sqf([1,-150,8500,-225000,2740000,-12000000], ZZ) == \
        [(QQ(10), QQ(10)), (QQ(20), QQ(20)), (QQ(30), QQ(30)), (QQ(40), QQ(40)), (QQ(50), QQ(50))]

    assert dup_isolate_real_roots_sqf([1,1], ZZ) == \
        [(-QQ(1), -QQ(1))]
    assert dup_isolate_real_roots_sqf([1,3,2], ZZ) == \
        [(-QQ(2), -QQ(2)), (-QQ(1), -QQ(1))]
    assert dup_isolate_real_roots_sqf([1,6,11,6], ZZ) == \
        [(-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1))]
    assert dup_isolate_real_roots_sqf([1,10,35,50,24], ZZ) == \
        [(-QQ(4), -QQ(4)), (-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1))]
    assert dup_isolate_real_roots_sqf([1,15,85,225,274,120], ZZ) == \
        [(-QQ(5), -QQ(5)), (-QQ(4), -QQ(4)), (-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1))]

    assert dup_isolate_real_roots_sqf([1,10], ZZ) == \
        [(-QQ(10), -QQ(10))]
    assert dup_isolate_real_roots_sqf([1,30,200], ZZ) == \
        [(-QQ(20), -QQ(20)), (-QQ(10), -QQ(10))]
    assert dup_isolate_real_roots_sqf([1,60,1100,6000], ZZ) == \
        [(-QQ(30), -QQ(30)), (-QQ(20), -QQ(20)), (-QQ(10), -QQ(10))]
    assert dup_isolate_real_roots_sqf([1,100,3500,50000,240000], ZZ) == \
        [(-QQ(40), -QQ(40)), (-QQ(30), -QQ(30)), (-QQ(20), -QQ(20)), (-QQ(10), -QQ(10))]
    assert dup_isolate_real_roots_sqf([1,150,8500,225000,2740000,12000000], ZZ) == \
        [(-QQ(50), -QQ(50)), (-QQ(40), -QQ(40)), (-QQ(30), -QQ(30)), (-QQ(20), -QQ(20)), (-QQ(10), -QQ(10))]

    assert dup_isolate_real_roots_sqf([1,0,-5], ZZ) == \
        [(QQ(-3), QQ(-2)), (QQ(2), QQ(3))]
    assert dup_isolate_real_roots_sqf([1,0,0,-5], ZZ) == \
        [(QQ(1), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,0,0,-5], ZZ) == \
        [(QQ(-2), QQ(-1)), (QQ(1), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,0,0,0,-5], ZZ) == \
        [(QQ(1), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,0,0,0,0,-5], ZZ) == \
        [(QQ(-2), QQ(-1)), (QQ(1), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,0,0,0,0,0,-5], ZZ) == \
        [(QQ(1), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,0,0,0,0,0,0,-5], ZZ) == \
        [(QQ(-2), QQ(-1)), (QQ(1), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,0,0,0,0,0,0,0,-5], ZZ) == \
        [(QQ(1), QQ(2))]

    assert dup_isolate_real_roots_sqf([1,0,-1], ZZ) == \
        [(-QQ(1), -QQ(1)), (QQ(1), QQ(1))]
    assert dup_isolate_real_roots_sqf([1,2,-1,-2], ZZ) == \
        [(-QQ(2), -QQ(2)), (-QQ(1), -QQ(1)), (QQ(1), QQ(1))]
    assert dup_isolate_real_roots_sqf([1,0,-5,0,4], ZZ) == \
        [(-QQ(2), -QQ(2)), (-QQ(1), -QQ(1)), (QQ(1), QQ(1)), (QQ(2), QQ(2))]
    assert dup_isolate_real_roots_sqf([1,3,-5,-15,4,12], ZZ) == \
        [(-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1)), (QQ(1), QQ(1)),
         ( QQ(2),  QQ(2))]
    assert dup_isolate_real_roots_sqf([1,0,-14,0,49,0,-36], ZZ) == \
        [(-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1)), (QQ(1), QQ(1)),
         ( QQ(2),  QQ(2)), ( QQ(3),  QQ(3))]
    assert dup_isolate_real_roots_sqf([2,1,-28,-14,98,49,-72,-36], ZZ) == \
        [(-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1)), (-QQ(1), QQ(0)),
         ( QQ(1),  QQ(1)), ( QQ(2),  QQ(2)), ( QQ(3),  QQ(3))]
    assert dup_isolate_real_roots_sqf([4,0,-57,0,210,0,-193,0,36], ZZ) == \
        [(-QQ(3), -QQ(3)), (-QQ(2), -QQ(2)), (-QQ(1), -QQ(1)), (-QQ(1), QQ(0)),
         ( QQ(0),  QQ(1)), ( QQ(1),  QQ(1)), ( QQ(2),  QQ(2)), ( QQ(3), QQ(3))]

    f = [9,0,-2]

    assert dup_isolate_real_roots_sqf(f, ZZ) == \
        [(QQ(-1), QQ(0)), (QQ(0), QQ(1))]

    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,10)) == \
        [(QQ(-1,2), QQ(-3,7)), (QQ(3,7), QQ(1,2))]
    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,100)) == \
        [(QQ(-9,19), QQ(-8,17)), (QQ(8,17), QQ(9,19))]
    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,1000)) == \
        [(QQ(-33,70), QQ(-8,17)), (QQ(8,17), QQ(33,70))]
    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,10000)) == \
        [(QQ(-33,70), QQ(-107,227)), (QQ(107,227), QQ(33,70))]
    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,100000)) == \
        [(QQ(-305,647), QQ(-272,577)), (QQ(272,577), QQ(305,647))]
    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,1000000)) == \
        [(QQ(-1121,2378), QQ(-272,577)), (QQ(272,577), QQ(1121,2378))]

    f = [200100012, -700390052, 700490079, -200240054, 40017, -2]

    assert dup_isolate_real_roots_sqf(f, ZZ) == \
        [(QQ(0), QQ(1,10002)), (QQ(1,10002), QQ(1,10002)), (QQ(1,2), QQ(1,2)), (QQ(1), QQ(1)), (QQ(2), QQ(2))]

    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,100000)) == \
        [(QQ(1,10003), QQ(1,10003)), (QQ(1,10002), QQ(1,10002)), (QQ(1,2), QQ(1,2)), (QQ(1), QQ(1)), (QQ(2), QQ(2))]

    a, b, c, d = 10000090000001, 2000100003, 10000300007, 10000005000008

    f = [ 20001600074001600021,
          1700135866278935491773999857,
         -2000179008931031182161141026995283662899200197,
         -800027600594323913802305066986600025,
          100000950000540000725000008]

    assert dup_isolate_real_roots_sqf(f, ZZ) == \
        [(-QQ(a), -QQ(a)), (-QQ(1,1),  QQ(0,1)), (QQ(0,1), QQ(1,1)), (QQ(d), QQ(d))]

    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,100000000000)) == \
        [(-QQ(a), -QQ(a)), (-QQ(1,b), -QQ(1,b)), (QQ(1,c), QQ(1,c)), (QQ(d), QQ(d))]

    (u, v), B, C, (s, t) = dup_isolate_real_roots_sqf(f, ZZ, fast=True)

    assert u < -a < v and B == (-QQ(1), QQ(0)) and C == (QQ(0), QQ(1)) and s < d < t

    assert dup_isolate_real_roots_sqf(f, ZZ, fast=True, eps=QQ(1,100000000000000000000000000000)) == \
        [(-QQ(a), -QQ(a)), (-QQ(1,b), -QQ(1,b)), (QQ(1,c), QQ(1,c)), (QQ(d), QQ(d))]

    assert dup_isolate_real_roots_sqf([QQ(8,5), QQ(-87374,3855), QQ(-17,771)], QQ) == \
        [(QQ(-1), QQ(0)), (QQ(14), QQ(15))]

    f = [-10, 8, 80, -32, -160]

    assert dup_isolate_real_roots_sqf(f, ZZ) == \
        [(-QQ(2), -QQ(2)), (-QQ(2), -QQ(1)), (QQ(2), QQ(2)), (QQ(2), QQ(3))]

    assert dup_isolate_real_roots_sqf(f, ZZ, eps=QQ(1,100)) == \
        [(-QQ(2), -QQ(2)), (-QQ(23,14), -QQ(18,11)), (QQ(2), QQ(2)), (QQ(39,16), QQ(22,9))]

    assert dup_isolate_real_roots_sqf([1, -1], ZZ, inf=2) == []
    assert dup_isolate_real_roots_sqf([1, -1], ZZ, sup=0) == []

    assert dup_isolate_real_roots_sqf([1, -1], ZZ) == [(1, 1)]
    assert dup_isolate_real_roots_sqf([1, -1], ZZ, inf=1) == [(1, 1)]
    assert dup_isolate_real_roots_sqf([1, -1], ZZ, sup=1) == [(1, 1)]
    assert dup_isolate_real_roots_sqf([1, -1], ZZ, inf=1, sup=1) == [(1, 1)]

    f = [1, 0, -2]

    assert dup_isolate_real_roots_sqf(f, ZZ, inf=QQ(7,4)) == []
    assert dup_isolate_real_roots_sqf(f, ZZ, inf=QQ(7,5)) == [(QQ(7,5), QQ(3,2))]
    assert dup_isolate_real_roots_sqf(f, ZZ, sup=QQ(7,5)) == [(-2, -1)]
    assert dup_isolate_real_roots_sqf(f, ZZ, sup=QQ(7,4)) == [(-2, -1), (1, QQ(3,2))]
    assert dup_isolate_real_roots_sqf(f, ZZ, sup=-QQ(7,4)) == []
    assert dup_isolate_real_roots_sqf(f, ZZ, sup=-QQ(7,5)) == [(-QQ(3,2), -QQ(7,5))]
    assert dup_isolate_real_roots_sqf(f, ZZ, inf=-QQ(7,5)) == [(1, 2)]
    assert dup_isolate_real_roots_sqf(f, ZZ, inf=-QQ(7,4)) == [(-QQ(3,2), -1), (1, 2)]

    I = [(-2, -1), (1, 2)]

    assert dup_isolate_real_roots_sqf(f, ZZ, inf=-2) == I
    assert dup_isolate_real_roots_sqf(f, ZZ, sup=+2) == I

    assert dup_isolate_real_roots_sqf(f, ZZ, inf=-2, sup=2) == I

    raises(DomainError, "dup_isolate_real_roots_sqf([EX(1), EX(2)], EX)")

def test_dup_isolate_real_roots():
    assert dup_isolate_real_roots([], ZZ) == []
    assert dup_isolate_real_roots([3], ZZ) == []

    assert dup_isolate_real_roots([5,0], ZZ) ==  [((QQ(0), QQ(0)), 1)]
    assert dup_isolate_real_roots([7,0,0,0,0], ZZ) == [((QQ(0), QQ(0)), 4)]

    assert dup_isolate_real_roots([1, 1,0], ZZ) == [((-QQ(1), -QQ(1)), 1), ((QQ(0), QQ(0)), 1)]
    assert dup_isolate_real_roots([1,-1,0], ZZ) == [(( QQ(0),  QQ(0)), 1), ((QQ(1), QQ(1)), 1)]

    assert dup_isolate_real_roots([1,0,0,1,1], ZZ) == []

    I = [((-QQ(2), -QQ(1)), 1), ((QQ(1), QQ(2)), 1)]

    assert dup_isolate_real_roots([1,0,-2], ZZ) == I
    assert dup_isolate_real_roots([-1,0,2], ZZ) == I

    f = [16,-96,24,936,-1599,-2880,9196,552,-21831,13968,21690,-26784,-2916,15552,-5832]
    g = dup_sqf_part(f, ZZ)

    assert dup_isolate_real_roots(f, ZZ) == \
        [((-QQ(2), -QQ(3,2)), 2), ((-QQ(3,2), -QQ(1,1)), 3),
         (( QQ(1),  QQ(3,2)), 3), (( QQ(3,2),  QQ(3,2)), 4), ((QQ(5,3), QQ(2)), 2)]

    assert dup_isolate_real_roots_sqf(g, ZZ) == \
        [(-QQ(2), -QQ(3,2)), (-QQ(3,2), -QQ(1,1)),
         ( QQ(1),  QQ(3,2)), ( QQ(3,2),  QQ(3,2)), (QQ(3,2), QQ(2))]
    assert dup_isolate_real_roots(g, ZZ) == \
        [((-QQ(2), -QQ(3,2)), 1), ((-QQ(3,2), -QQ(1,1)), 1),
         (( QQ(1),  QQ(3,2)), 1), (( QQ(3,2),  QQ(3,2)), 1), ((QQ(3,2), QQ(2)), 1)]

    assert dup_isolate_real_roots([1, -1], ZZ, inf=2) == []
    assert dup_isolate_real_roots([1, -1], ZZ, sup=0) == []

    assert dup_isolate_real_roots([1, -1], ZZ) == [((1, 1), 1)]
    assert dup_isolate_real_roots([1, -1], ZZ, inf=1) == [((1, 1), 1)]
    assert dup_isolate_real_roots([1, -1], ZZ, sup=1) == [((1, 1), 1)]
    assert dup_isolate_real_roots([1, -1], ZZ, inf=1, sup=1) == [((1, 1), 1)]

    f = [1, 0, -4, 0, 4]

    assert dup_isolate_real_roots(f, ZZ, inf=QQ(7,4)) == []
    assert dup_isolate_real_roots(f, ZZ, inf=QQ(7,5)) == [((QQ(7,5), QQ(3,2)), 2)]
    assert dup_isolate_real_roots(f, ZZ, sup=QQ(7,5)) == [((-2, -1), 2)]
    assert dup_isolate_real_roots(f, ZZ, sup=QQ(7,4)) == [((-2, -1), 2), ((1, QQ(3,2)), 2)]
    assert dup_isolate_real_roots(f, ZZ, sup=-QQ(7,4)) == []
    assert dup_isolate_real_roots(f, ZZ, sup=-QQ(7,5)) == [((-QQ(3,2), -QQ(7,5)), 2)]
    assert dup_isolate_real_roots(f, ZZ, inf=-QQ(7,5)) == [((1, 2), 2)]
    assert dup_isolate_real_roots(f, ZZ, inf=-QQ(7,4)) == [((-QQ(3,2), -1), 2), ((1, 2), 2)]

    I = [((-2, -1), 2), ((1, 2), 2)]

    assert dup_isolate_real_roots(f, ZZ, inf=-2) == I
    assert dup_isolate_real_roots(f, ZZ, sup=+2) == I

    assert dup_isolate_real_roots(f, ZZ, inf=-2, sup=2) == I

    f = [1, -3, -1, 11, -8, -8, 12, -4, 0, 0, 0, 0]

    assert dup_isolate_real_roots(f, ZZ, basis=False) == \
        [((-2, -1), 2), ((0, 0), 4), ((1, 1), 3), ((1, 2), 2)]
    assert dup_isolate_real_roots(f, ZZ, basis=True) == \
        [((-2, -1), 2, [1, 0, -2]), ((0, 0), 4, [1, 0]), ((1, 1), 3, [1, -1]), ((1, 2), 2, [1, 0, -2])]

    raises(DomainError, "dup_isolate_real_roots([EX(1), EX(2)], EX)")

def test_dup_isolate_real_roots_list():
    assert dup_isolate_real_roots_list([[1, 1,0], [1,0]], ZZ) == \
        [((-QQ(1), -QQ(1)), {0: 1}), ((QQ(0), QQ(0)), {0: 1, 1: 1})]
    assert dup_isolate_real_roots_list([[1,-1,0], [1,0]], ZZ) == \
        [((QQ(0), QQ(0)), {0: 1, 1: 1}), ((QQ(1), QQ(1)), {0: 1})]

    f = dup_from_raw_dict({5: ZZ(1), 0: -ZZ(200)}, ZZ)
    g = dup_from_raw_dict({5: ZZ(1), 0: -ZZ(201)}, ZZ)

    assert dup_isolate_real_roots_list([f, g], ZZ) == \
        [((QQ(75,26), QQ(101,35)), {0: 1}), ((QQ(283,98), QQ(26,9)), {1: 1})]

    f = dup_from_raw_dict({5: -QQ(1,200), 0: QQ(1)}, QQ)
    g = dup_from_raw_dict({5: -QQ(1,201), 0: QQ(1)}, QQ)

    assert dup_isolate_real_roots_list([f, g], QQ) == \
        [((QQ(75,26), QQ(101,35)), {0: 1}), ((QQ(283,98), QQ(26,9)), {1: 1})]

    assert dup_isolate_real_roots_list([[1,1], [1,2], [1,-1], [1,1], [1,-1], [1,-1]], ZZ) == \
        [((-QQ(2), -QQ(2)), {1: 1}), ((-QQ(1), -QQ(1)), {0: 1, 3: 1}), ((QQ(1), QQ(1)), {2: 1, 4: 1, 5: 1})]

    assert dup_isolate_real_roots_list([[1,1], [1,2], [1,-1], [1,1], [1,-1], [1,2]], ZZ) == \
        [((-QQ(2), -QQ(2)), {1: 1, 5: 1}), ((-QQ(1), -QQ(1)), {0: 1, 3: 1}), ((QQ(1), QQ(1)), {2: 1, 4: 1})]

    f, g = [1, 0, -4, 0, 4], [1, -1]

    assert dup_isolate_real_roots_list([f, g], ZZ, inf=QQ(7,4)) == []
    assert dup_isolate_real_roots_list([f, g], ZZ, inf=QQ(7,5)) == [((QQ(7,5), QQ(3,2)), {0: 2})]
    assert dup_isolate_real_roots_list([f, g], ZZ, sup=QQ(7,5)) == [((-2, -1), {0: 2}), ((1, 1), {1: 1})]
    assert dup_isolate_real_roots_list([f, g], ZZ, sup=QQ(7,4)) == [((-2, -1), {0: 2}), ((1, 1), {1: 1}), ((1, QQ(3,2)), {0: 2})]
    assert dup_isolate_real_roots_list([f, g], ZZ, sup=-QQ(7,4)) == []
    assert dup_isolate_real_roots_list([f, g], ZZ, sup=-QQ(7,5)) == [((-QQ(3,2), -QQ(7,5)), {0: 2})]
    assert dup_isolate_real_roots_list([f, g], ZZ, inf=-QQ(7,5)) == [((1, 1), {1: 1}), ((1, 2), {0: 2})]
    assert dup_isolate_real_roots_list([f, g], ZZ, inf=-QQ(7,4)) == [((-QQ(3,2), -1), {0: 2}), ((1, 1), {1: 1}), ((1, 2), {0: 2})]

    f, g = [2, 0, -1], [1, 0, -2]

    assert dup_isolate_real_roots_list([f, g], ZZ) == \
        [((-QQ(2), -QQ(1)), {1: 1}), ((-QQ(1), QQ(0)), {0: 1}), ((QQ(0), QQ(1)), {0: 1}), ((QQ(1), QQ(2)), {1: 1})]
    assert dup_isolate_real_roots_list([f, g], ZZ, strict=True) == \
        [((-QQ(3,2), -QQ(4,3)), {1: 1}), ((-QQ(1), -QQ(2,3)), {0: 1}), ((QQ(2,3), QQ(1)), {0: 1}), ((QQ(4,3), QQ(3,2)), {1: 1})]

    f, g = [1, 0, -2], [1, -1, -2, 2]

    assert dup_isolate_real_roots_list([f, g], ZZ) == \
        [((-QQ(2), -QQ(1)), {1: 1, 0: 1}), ((QQ(1), QQ(1)), {1: 1}), ((QQ(1), QQ(2)), {1: 1, 0: 1})]

    f, g = [1, 0, -2, 0], [1, -1, -2, 2, 0, 0]

    assert dup_isolate_real_roots_list([f, g], ZZ) == \
        [((-QQ(2), -QQ(1)), {1: 1, 0: 1}), ((QQ(0), QQ(0)), {0: 1, 1: 2}), ((QQ(1), QQ(1)), {1: 1}), ((QQ(1), QQ(2)), {1: 1, 0: 1})]

    f, g = [1, -3, -1, 11, -8, -8, 12, -4, 0, 0], [1, -2, 3, -4, 2, 0]

    assert dup_isolate_real_roots_list([f, g], ZZ, basis=False) == \
        [((-2, -1), {0: 2}), ((0, 0), {0: 2, 1: 1}), ((1, 1), {0: 3, 1: 2}), ((1, 2), {0: 2})]
    assert dup_isolate_real_roots_list([f, g], ZZ, basis=True) == \
        [((-2, -1), {0: 2}, [1, 0, -2]), ((0, 0), {0: 2, 1: 1}, [1, 0]), ((1, 1), {0: 3, 1: 2}, [1, -1]), ((1, 2), {0: 2}, [1, 0, -2])]

    raises(DomainError, "dup_isolate_real_roots_list([[EX(1), EX(2)]], EX)")

def test_dup_count_real_roots():
    assert dup_count_real_roots([], ZZ) == 0
    assert dup_count_real_roots([7], ZZ) == 0

    assert dup_count_real_roots([1,-1], ZZ) == 1
    assert dup_count_real_roots([1,-1], ZZ, inf=1) == 1
    assert dup_count_real_roots([1,-1], ZZ, sup=0) == 0
    assert dup_count_real_roots([1,-1], ZZ, sup=1) == 1
    assert dup_count_real_roots([1,-1], ZZ, inf=0, sup=1) == 1
    assert dup_count_real_roots([1,-1], ZZ, inf=0, sup=2) == 1
    assert dup_count_real_roots([1,-1], ZZ, inf=1, sup=2) == 1

    assert dup_count_real_roots([1,0,-2], ZZ) == 2
    assert dup_count_real_roots([1,0,-2], ZZ, sup=0) == 1
    assert dup_count_real_roots([1,0,-2], ZZ, inf=-1, sup=1) == 0

a, b = (-QQ(1), -QQ(1)), (QQ(1), QQ(1))
c, d = ( QQ(0),  QQ(0)), (QQ(1), QQ(1))

def test_dup_count_complex_roots_1():
    # z-1
    assert dup_count_complex_roots([1,-1], ZZ, a, b) == 1
    assert dup_count_complex_roots([1,-1], ZZ, c, d) == 1

    # z+1
    assert dup_count_complex_roots([1,1], ZZ, a, b) == 1
    assert dup_count_complex_roots([1,1], ZZ, c, d) == 0

def test_dup_count_complex_roots_2():
    # (z-1)*(z)
    assert dup_count_complex_roots([1,-1,0], ZZ, a, b) == 2
    assert dup_count_complex_roots([1,-1,0], ZZ, c, d) == 2

    # (z-1)*(-z)
    assert dup_count_complex_roots([-1,1,0], ZZ, a, b) == 2
    assert dup_count_complex_roots([-1,1,0], ZZ, c, d) == 2

    # (z+1)*(z)
    assert dup_count_complex_roots([1,1,0], ZZ, a, b) == 2
    assert dup_count_complex_roots([1,1,0], ZZ, c, d) == 1

    # (z+1)*(-z)
    assert dup_count_complex_roots([-1,-1,0], ZZ, a, b) == 2
    assert dup_count_complex_roots([-1,-1,0], ZZ, c, d) == 1

def test_dup_count_complex_roots_3():
    # (z-1)*(z+1)
    assert dup_count_complex_roots([1,0,-1], ZZ, a, b) == 2
    assert dup_count_complex_roots([1,0,-1], ZZ, c, d) == 1

    # (z-1)*(z+1)*(z)
    assert dup_count_complex_roots([1,0,-1,0], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,0,-1,0], ZZ, c, d) == 2

    # (z-1)*(z+1)*(-z)
    assert dup_count_complex_roots([-1,0,1,0], ZZ, a, b) == 3
    assert dup_count_complex_roots([-1,0,1,0], ZZ, c, d) == 2

def test_dup_count_complex_roots_4():
    # (z-I)*(z+I)
    assert dup_count_complex_roots([1,0,1], ZZ, a, b) == 2
    assert dup_count_complex_roots([1,0,1], ZZ, c, d) == 1

    # (z-I)*(z+I)*(z)
    assert dup_count_complex_roots([1,0,1,0], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,0,1,0], ZZ, c, d) == 2

    # (z-I)*(z+I)*(-z)
    assert dup_count_complex_roots([-1,0,-1,0], ZZ, a, b) == 3
    assert dup_count_complex_roots([-1,0,-1,0], ZZ, c, d) == 2

    # (z-I)*(z+I)*(z-1)
    assert dup_count_complex_roots([1,-1,1,-1], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,-1,1,-1], ZZ, c, d) == 2

    # (z-I)*(z+I)*(z-1)*(z)
    assert dup_count_complex_roots([1,-1,1,-1,0], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,-1,1,-1,0], ZZ, c, d) == 3

    # (z-I)*(z+I)*(z-1)*(-z)
    assert dup_count_complex_roots([-1,1,-1,1,0], ZZ, a, b) == 4
    assert dup_count_complex_roots([-1,1,-1,1,0], ZZ, c, d) == 3

    # (z-I)*(z+I)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,0,0,0,-1], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,0,0,0,-1], ZZ, c, d) == 2

    # (z-I)*(z+I)*(z-1)*(z+1)*(z)
    assert dup_count_complex_roots([1,0,0,0,-1,0], ZZ, a, b) == 5
    assert dup_count_complex_roots([1,0,0,0,-1,0], ZZ, c, d) == 3

    # (z-I)*(z+I)*(z-1)*(z+1)*(-z)
    assert dup_count_complex_roots([-1,0,0,0,1,0], ZZ, a, b) == 5
    assert dup_count_complex_roots([-1,0,0,0,1,0], ZZ, c, d) == 3

def test_dup_count_complex_roots_5():
    # (z-I+1)*(z+I+1)
    assert dup_count_complex_roots([1,2,2], ZZ, a, b) == 2
    assert dup_count_complex_roots([1,2,2], ZZ, c, d) == 0

    # (z-I+1)*(z+I+1)*(z-1)
    assert dup_count_complex_roots([1,1,0,-2], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,1,0,-2], ZZ, c, d) == 1

    # (z-I+1)*(z+I+1)*(z-1)*z
    assert dup_count_complex_roots([1,1,0,-2,0], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,1,0,-2,0], ZZ, c, d) == 2

    # (z-I+1)*(z+I+1)*(z+1)
    assert dup_count_complex_roots([1,3,4,2], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,3,4,2], ZZ, c, d) == 0

    # (z-I+1)*(z+I+1)*(z+1)*z
    assert dup_count_complex_roots([1,3,4,2,0], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,3,4,2,0], ZZ, c, d) == 1

    # (z-I+1)*(z+I+1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,2,1,-2,-2], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,2,1,-2,-2], ZZ, c, d) == 1

    # (z-I+1)*(z+I+1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,2,1,-2,-2,0], ZZ, a, b) == 5
    assert dup_count_complex_roots([1,2,1,-2,-2,0], ZZ, c, d) == 2

def test_dup_count_complex_roots_6():
    # (z-I-1)*(z+I-1)
    assert dup_count_complex_roots([1,-2,2], ZZ, a, b) == 2
    assert dup_count_complex_roots([1,-2,2], ZZ, c, d) == 1

    # (z-I-1)*(z+I-1)*(z-1)
    assert dup_count_complex_roots([1,-3,4,-2], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,-3,4,-2], ZZ, c, d) == 2

    # (z-I-1)*(z+I-1)*(z-1)*z
    assert dup_count_complex_roots([1,-3,4,-2,0], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,-3,4,-2,0], ZZ, c, d) == 3

    # (z-I-1)*(z+I-1)*(z+1)
    assert dup_count_complex_roots([1,-1,0,2], ZZ, a, b) == 3
    assert dup_count_complex_roots([1,-1,0,2], ZZ, c, d) == 1

    # (z-I-1)*(z+I-1)*(z+1)*z
    assert dup_count_complex_roots([1,-1,0,2,0], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,-1,0,2,0], ZZ, c, d) == 2

    # (z-I-1)*(z+I-1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,-2,1,2,-2], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,-2,1,2,-2], ZZ, c, d) == 2

    # (z-I-1)*(z+I-1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,-2,1,2,-2,0], ZZ, a, b) == 5
    assert dup_count_complex_roots([1,-2,1,2,-2,0], ZZ, c, d) == 3

def test_dup_count_complex_roots_7():
    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)
    assert dup_count_complex_roots([1,0,0,0,4], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,0,0,0,4], ZZ, c, d) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-2)
    assert dup_count_complex_roots([1,-2,0,0,4,-8], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,-2,0,0,4,-8], ZZ, c, d) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z**2-2)
    assert dup_count_complex_roots([1,0,-2,0,4,0,-8], ZZ, a, b) == 4
    assert dup_count_complex_roots([1,0,-2,0,4,0,-8], ZZ, c, d) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)
    assert dup_count_complex_roots([1,-1,0,0,4,-4], ZZ, a, b) == 5
    assert dup_count_complex_roots([1,-1,0,0,4,-4], ZZ, c, d) == 2

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*z
    assert dup_count_complex_roots([1,-1,0,0,4,-4,0], ZZ, a, b) == 6
    assert dup_count_complex_roots([1,-1,0,0,4,-4,0], ZZ, c, d) == 3

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z+1)
    assert dup_count_complex_roots([1,1,0,0,4,4], ZZ, a, b) == 5
    assert dup_count_complex_roots([1,1,0,0,4,4], ZZ, c, d) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z+1)*z
    assert dup_count_complex_roots([1,1,0,0,4,4,0], ZZ, a, b) == 6
    assert dup_count_complex_roots([1,1,0,0,4,4,0], ZZ, c, d) == 2

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4], ZZ, a, b) == 6
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4], ZZ, c, d) == 2

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4,0], ZZ, a, b) == 7
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4,0], ZZ, c, d) == 3

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4], ZZ, a, b) == 8
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4], ZZ, c, d) == 3

def test_dup_count_complex_roots_8():
    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)*z
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4,0], ZZ, a, b) == 9
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4,0], ZZ, c, d) == 4

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)*(z**2-2)*z
    assert dup_count_complex_roots([1,0,-2,0,3,0,-6,0,-4,0,8,0], ZZ, a, b) == 9
    assert dup_count_complex_roots([1,0,-2,0,3,0,-6,0,-4,0,8,0], ZZ, c, d) == 4

def test_dup_count_complex_roots_implicit():
    f = [1,0,0,0,-1,0] # z*(z-1)*(z+1)*(z-I)*(z+I)

    assert dup_count_complex_roots(f, ZZ) == 5

    assert dup_count_complex_roots(f, ZZ, sup=(0, 0)) == 3
    assert dup_count_complex_roots(f, ZZ, inf=(0, 0)) == 3

def test_dup_count_complex_roots_exclude():
    f = [1,0,0,0,-1,0] # z*(z-1)*(z+1)*(z-I)*(z+I)

    a, b = (-QQ(1), QQ(0)), (QQ(1), QQ(1))

    assert dup_count_complex_roots(f, ZZ, a, b) == 4

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['S']) == 3
    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['N']) == 3

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['S', 'N']) == 2

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['E']) == 4
    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['W']) == 4

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['E', 'W']) == 4

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['N', 'S', 'E', 'W']) == 2

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['SW']) == 3
    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['SE']) == 3

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['SW', 'SE']) == 2
    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['SW', 'SE', 'S']) == 1
    assert dup_count_complex_roots(f, ZZ, a, b, exclude=['SW', 'SE', 'S', 'N']) == 0

    a, b = (QQ(0), QQ(0)), (QQ(1), QQ(1))

    assert dup_count_complex_roots(f, ZZ, a, b, exclude=True) == 1

def test_dup_isolate_complex_roots_sqf():
    f = [1, -2, 3]

    assert dup_isolate_complex_roots_sqf(f, ZZ) == \
        [((0, 0), (6, 6)), ((0, -6), (6, 0))]

    assert dup_isolate_complex_roots_sqf(f, ZZ, eps=QQ(1,10)) == \
        [((QQ(15,16), QQ(45,32)), (QQ(33,32), QQ(3,2))), ((QQ(15,16), -QQ(3,2)), (QQ(33,32), -QQ(45,32)))]
    assert dup_isolate_complex_roots_sqf(f, ZZ, eps=QQ(1,100)) == \
        [((QQ(255,256), QQ(723,512)), (QQ(513,512), QQ(363,256))), ((QQ(255,256), -QQ(363,256)), (QQ(513,512), -QQ(723,512)))]

    f = [7, -19, 20, 17, 20]

    assert dup_isolate_complex_roots_sqf(f, ZZ) == \
        [((-QQ(40,7), 0), (0, QQ(40,7))), ((-QQ(40,7), -QQ(40,7)), (0, 0)), ((0, 0), (QQ(40,7), QQ(40,7))), ((0, -QQ(40,7)), (QQ(40,7), 0))]

def test_dup_isolate_all_roots_sqf():
    f = [4, -1, 2, 5, 0]

    assert dup_isolate_all_roots_sqf(f, ZZ) == \
        ([(-1, 0), (0, 0)], [((0, 0), (QQ(5,2), QQ(5,2))), ((0, -QQ(5,2)), (QQ(5,2), 0))])

    assert dup_isolate_all_roots_sqf(f, ZZ, eps=QQ(1,10)) == \
        ([(QQ(-7,8), QQ(-6,7)), (0, 0)], [((QQ(35,64), QQ(65,64)), (QQ(5,8), QQ(35,32))), ((QQ(35,64), -QQ(35,32)), (QQ(5,8), -QQ(65,64)))])

def test_dup_isolate_all_roots():
    f = [4, -1, 2, 5, 0]

    assert dup_isolate_all_roots(f, ZZ) == \
        ([((-1, 0), 1), ((0, 0), 1)], [(((0, 0), (QQ(5,2), QQ(5,2))), 1), (((0, -QQ(5,2)), (QQ(5,2), 0)), 1)])

    assert dup_isolate_all_roots(f, ZZ, eps=QQ(1,10)) == \
        ([((QQ(-7,8), QQ(-6,7)), 1), ((0, 0), 1)], [(((QQ(35,64), QQ(65,64)), (QQ(5,8), QQ(35,32))), 1), (((QQ(35,64), -QQ(35,32)), (QQ(5,8), -QQ(65,64))), 1)])

    raises(NotImplementedError, "dup_isolate_all_roots([1, 1, -2, -2, 1, 1], ZZ)")

