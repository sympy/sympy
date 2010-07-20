"""Tests for dense recursive polynomials' arithmetics. """

from sympy.polys.densebasic import (
    dup_normal, dmp_normal,
)

from sympy.polys.densearith import (
    dup_add_term, dmp_add_term,
    dup_sub_term, dmp_sub_term,
    dup_mul_term, dmp_mul_term,
    dup_mul_ground, dmp_mul_ground,
    dup_quo_ground, dmp_quo_ground,
    dup_exquo_ground, dmp_exquo_ground,
    dup_lshift, dup_rshift,
    dup_abs, dmp_abs,
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_sqr, dmp_sqr,
    dup_pow, dmp_pow,
    dup_add_mul, dmp_add_mul,
    dup_sub_mul, dmp_sub_mul,
    dup_pdiv, dup_prem, dup_pquo, dup_pexquo,
    dmp_pdiv, dmp_prem, dmp_pquo, dmp_pexquo,
    dup_rr_div, dmp_rr_div,
    dup_ff_div, dmp_ff_div,
    dup_div, dup_rem, dup_quo, dup_exquo,
    dmp_div, dmp_rem, dmp_quo, dmp_exquo,
    dup_max_norm, dmp_max_norm,
    dup_l1_norm, dmp_l1_norm,
    dup_expand, dmp_expand,
    dup_revert,
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed, NotReversible,
)

from sympy.polys.specialpolys import f_0
from sympy.polys.domains import FF, ZZ, QQ

from sympy import raises

F_0 = dmp_mul_ground(dmp_normal(f_0, 2, QQ), QQ(1,7), 2, QQ)

def test_dup_add_term():
    f = dup_normal([], ZZ)

    assert dup_add_term(f, ZZ(0), 0, ZZ) == dup_normal([], ZZ)

    assert dup_add_term(f, ZZ(1), 0, ZZ) == dup_normal([1], ZZ)
    assert dup_add_term(f, ZZ(1), 1, ZZ) == dup_normal([1, 0], ZZ)
    assert dup_add_term(f, ZZ(1), 2, ZZ) == dup_normal([1, 0, 0], ZZ)

    f = dup_normal([1,1,1], ZZ)

    assert dup_add_term(f, ZZ(1), 0, ZZ) == dup_normal([1, 1, 2], ZZ)
    assert dup_add_term(f, ZZ(1), 1, ZZ) == dup_normal([1, 2, 1], ZZ)
    assert dup_add_term(f, ZZ(1), 2, ZZ) == dup_normal([2, 1, 1], ZZ)

    assert dup_add_term(f, ZZ(1), 3, ZZ) == dup_normal([1, 1, 1, 1], ZZ)
    assert dup_add_term(f, ZZ(1), 4, ZZ) == dup_normal([1, 0, 1, 1, 1], ZZ)
    assert dup_add_term(f, ZZ(1), 5, ZZ) == dup_normal([1, 0, 0, 1, 1, 1], ZZ)
    assert dup_add_term(f, ZZ(1), 6, ZZ) == dup_normal([1, 0, 0, 0, 1, 1, 1], ZZ)

    assert dup_add_term(f,ZZ(-1), 2, ZZ) == dup_normal([1, 1], ZZ)

def test_dmp_add_term():
    assert dmp_add_term([ZZ(1),ZZ(1),ZZ(1)], ZZ(1), 2, 0, ZZ) == \
           dup_add_term([ZZ(1),ZZ(1),ZZ(1)], ZZ(1), 2, ZZ)
    assert dmp_add_term(f_0, [[]], 3, 2, ZZ) == f_0
    assert dmp_add_term(F_0, [[]], 3, 2, QQ) == F_0

def test_dup_sub_term():
    f = dup_normal([], ZZ)

    assert dup_sub_term(f, ZZ(0), 0, ZZ) == dup_normal([], ZZ)

    assert dup_sub_term(f, ZZ(1), 0, ZZ) == dup_normal([-1], ZZ)
    assert dup_sub_term(f, ZZ(1), 1, ZZ) == dup_normal([-1, 0], ZZ)
    assert dup_sub_term(f, ZZ(1), 2, ZZ) == dup_normal([-1, 0, 0], ZZ)

    f = dup_normal([1,1,1], ZZ)

    assert dup_sub_term(f, ZZ(2), 0, ZZ) == dup_normal([ 1, 1,-1], ZZ)
    assert dup_sub_term(f, ZZ(2), 1, ZZ) == dup_normal([ 1,-1, 1], ZZ)
    assert dup_sub_term(f, ZZ(2), 2, ZZ) == dup_normal([-1, 1, 1], ZZ)

    assert dup_sub_term(f, ZZ(1), 3, ZZ) == dup_normal([-1, 1, 1, 1], ZZ)
    assert dup_sub_term(f, ZZ(1), 4, ZZ) == dup_normal([-1, 0, 1, 1, 1], ZZ)
    assert dup_sub_term(f, ZZ(1), 5, ZZ) == dup_normal([-1, 0, 0, 1, 1, 1], ZZ)
    assert dup_sub_term(f, ZZ(1), 6, ZZ) == dup_normal([-1, 0, 0, 0, 1, 1, 1], ZZ)

    assert dup_sub_term(f, ZZ(1), 2, ZZ) == dup_normal([1, 1], ZZ)

def test_dmp_sub_term():
    assert dmp_sub_term([ZZ(1),ZZ(1),ZZ(1)], ZZ(1), 2, 0, ZZ) == \
           dup_sub_term([ZZ(1),ZZ(1),ZZ(1)], ZZ(1), 2, ZZ)
    assert dmp_sub_term(f_0, [[]], 3, 2, ZZ) == f_0
    assert dmp_sub_term(F_0, [[]], 3, 2, QQ) == F_0

def test_dup_mul_term():
    f = dup_normal([], ZZ)

    assert dup_mul_term(f, ZZ(2), 3, ZZ) == dup_normal([], ZZ)

    f = dup_normal([1,1], ZZ)

    assert dup_mul_term(f, ZZ(0), 3, ZZ) == dup_normal([], ZZ)

    f = dup_normal([1,2,3], ZZ)

    assert dup_mul_term(f, ZZ(2), 0, ZZ) == dup_normal([2,4,6], ZZ)
    assert dup_mul_term(f, ZZ(2), 1, ZZ) == dup_normal([2,4,6,0], ZZ)
    assert dup_mul_term(f, ZZ(2), 2, ZZ) == dup_normal([2,4,6,0,0], ZZ)
    assert dup_mul_term(f, ZZ(2), 3, ZZ) == dup_normal([2,4,6,0,0,0], ZZ)

def test_dmp_mul_term():
    assert dmp_mul_term([ZZ(1),ZZ(2),ZZ(3)], ZZ(2), 1, 0, ZZ) == \
           dup_mul_term([ZZ(1),ZZ(2),ZZ(3)], ZZ(2), 1, ZZ)

    assert dmp_mul_term([[]], [ZZ(2)], 3, 1, ZZ) == [[]]
    assert dmp_mul_term([[ZZ(1)]], [], 3, 1, ZZ) == [[]]

    assert dmp_mul_term([[ZZ(1),ZZ(2)], [ZZ(3)]], [ZZ(2)], 2, 1, ZZ) == \
               [[ZZ(2),ZZ(4)], [ZZ(6)], [], []]

    assert dmp_mul_term([[]], [QQ(2,3)], 3, 1, QQ) == [[]]
    assert dmp_mul_term([[QQ(1,2)]], [], 3, 1, QQ) == [[]]

    assert dmp_mul_term([[QQ(1,5),QQ(2,5)], [QQ(3,5)]], [QQ(2,3)], 2, 1, QQ) == \
               [[QQ(2,15),QQ(4,15)], [QQ(6,15)], [], []]

def test_dup_mul_ground():
    f = dup_normal([], ZZ)

    assert dup_mul_ground(f, ZZ(2), ZZ) == dup_normal([], ZZ)

    f = dup_normal([1,2,3], ZZ)

    assert dup_mul_ground(f, ZZ(0), ZZ) == dup_normal([], ZZ)
    assert dup_mul_ground(f, ZZ(2), ZZ) == dup_normal([2,4,6], ZZ)

def test_dmp_mul_ground():
    assert dmp_mul_ground(f_0, ZZ(2), 2, ZZ) == [
        [[ZZ(2),ZZ(4),ZZ(6)], [ZZ(4)]],
        [[ZZ(6)]],
        [[ZZ(8),ZZ(10),ZZ(12)], [ZZ(2),ZZ(4),ZZ(2)], [ZZ(2)]]
    ]

    assert dmp_mul_ground(F_0, QQ(1,2), 2, QQ) == [
        [[QQ(1,14),QQ(2,14),QQ(3,14)], [QQ(2,14)]],
        [[QQ(3,14)]],
        [[QQ(4,14),QQ(5,14),QQ(6,14)], [QQ(1,14),QQ(2,14),QQ(1,14)], [QQ(1,14)]]
    ]

def test_dup_quo_ground():
    raises(ZeroDivisionError, 'dup_quo_ground(dup_normal([1,2,3], ZZ), ZZ(0), ZZ)')
    raises(ExactQuotientFailed, 'dup_quo_ground(dup_normal([1,2,3], ZZ), ZZ(3), ZZ)')

    f = dup_normal([], ZZ)

    assert dup_quo_ground(f, ZZ(3), ZZ) == dup_normal([], ZZ)

    f = dup_normal([6,2,8], ZZ)

    assert dup_quo_ground(f, ZZ(1), ZZ) == f
    assert dup_quo_ground(f, ZZ(2), ZZ) == dup_normal([3,1,4], ZZ)

    f = dup_normal([6,2,8], QQ)

    assert dup_quo_ground(f, QQ(1), QQ) == f
    assert dup_quo_ground(f, QQ(2), QQ) == [QQ(3),QQ(1),QQ(4)]
    assert dup_quo_ground(f, QQ(7), QQ) == [QQ(6,7),QQ(2,7),QQ(8,7)]

def test_dup_exquo_ground():
    raises(ZeroDivisionError, 'dup_exquo_ground(dup_normal([1,2,3], ZZ), ZZ(0), ZZ)')

    f = dup_normal([], ZZ)

    assert dup_quo_ground(f, ZZ(3), ZZ) == dup_normal([], ZZ)

    f = dup_normal([6,2,8], ZZ)

    assert dup_exquo_ground(f, ZZ(1), ZZ) == f
    assert dup_exquo_ground(f, ZZ(2), ZZ) == dup_normal([3,1,4], ZZ)

    assert dup_exquo_ground(f, ZZ(3), ZZ) == dup_normal([2,0,2], ZZ)

    f = dup_normal([6,2,8], QQ)

    assert dup_exquo_ground(f, QQ(1), QQ) == f
    assert dup_exquo_ground(f, QQ(2), QQ) == [QQ(3),QQ(1),QQ(4)]
    assert dup_exquo_ground(f, QQ(7), QQ) == [QQ(6,7),QQ(2,7),QQ(8,7)]

def test_dmp_quo_ground():
    f = dmp_normal([[6],[2],[8]], 1, ZZ)

    assert dmp_quo_ground(f, ZZ(1), 1, ZZ) == f
    assert dmp_quo_ground(f, ZZ(2), 1, ZZ) == dmp_normal([[3],[1],[4]], 1, ZZ)

def test_dmp_exquo_ground():
    f = dmp_normal([[6],[2],[8]], 1, ZZ)

    assert dmp_exquo_ground(f, ZZ(1), 1, ZZ) == f
    assert dmp_exquo_ground(f, ZZ(2), 1, ZZ) == dmp_normal([[3],[1],[4]], 1, ZZ)

    assert dmp_normal(dmp_exquo_ground(f, ZZ(3), 1, ZZ), 1, ZZ) == dmp_normal([[2],[],[2]], 1, ZZ)

def test_dup_lshift():
    assert dup_lshift([], 3, ZZ) == []
    assert dup_lshift([1], 3, ZZ) == [1,0,0,0]

def test_dup_rshift():
    assert dup_rshift([], 3, ZZ) == []
    assert dup_rshift([1,0,0,0], 3, ZZ) == [1]

def test_dup_abs():
    assert dup_abs([], ZZ) == []
    assert dup_abs([ZZ( 1)], ZZ) == [ZZ(1)]
    assert dup_abs([ZZ(-7)], ZZ) == [ZZ(7)]
    assert dup_abs([ZZ(-1),ZZ(2),ZZ(3)], ZZ) == [ZZ(1),ZZ(2),ZZ(3)]

    assert dup_abs([], QQ) == []
    assert dup_abs([QQ( 1,2)], QQ) == [QQ(1,2)]
    assert dup_abs([QQ(-7,3)], QQ) == [QQ(7,3)]
    assert dup_abs([QQ(-1,7),QQ(2,7),QQ(3,7)], QQ) == [QQ(1,7),QQ(2,7),QQ(3,7)]

def test_dmp_abs():
    assert dmp_abs([ZZ(-1)], 0, ZZ) == [ZZ(1)]
    assert dmp_abs([QQ(-1,2)], 0, QQ) == [QQ(1,2)]

    assert dmp_abs([[[]]], 2, ZZ) == [[[]]]
    assert dmp_abs([[[ZZ(1)]]], 2, ZZ) == [[[ZZ(1)]]]
    assert dmp_abs([[[ZZ(-7)]]], 2, ZZ) == [[[ZZ(7)]]]

    assert dmp_abs([[[]]], 2, QQ) == [[[]]]
    assert dmp_abs([[[QQ(1,2)]]], 2, QQ) == [[[QQ(1,2)]]]
    assert dmp_abs([[[QQ(-7,9)]]], 2, QQ) == [[[QQ(7,9)]]]

def test_dup_neg():
    assert dup_neg([], ZZ) == []
    assert dup_neg([ZZ(1)], ZZ) == [ZZ(-1)]
    assert dup_neg([ZZ(-7)], ZZ) == [ZZ(7)]
    assert dup_neg([ZZ(-1),ZZ(2),ZZ(3)], ZZ) == [ZZ(1),ZZ(-2),ZZ(-3)]

    assert dup_neg([], QQ) == []
    assert dup_neg([QQ(1,2)], QQ) == [QQ(-1,2)]
    assert dup_neg([QQ(-7,9)], QQ) == [QQ(7,9)]
    assert dup_neg([QQ(-1,7),QQ(2,7),QQ(3,7)], QQ) == [QQ(1,7),QQ(-2,7),QQ(-3,7)]

def test_dmp_neg():
    assert dmp_neg([ZZ(-1)], 0, ZZ) == [ZZ(1)]
    assert dmp_neg([QQ(-1,2)], 0, QQ) == [QQ(1,2)]

    assert dmp_neg([[[]]], 2, ZZ) == [[[]]]
    assert dmp_neg([[[ZZ(1)]]], 2, ZZ) == [[[ZZ(-1)]]]
    assert dmp_neg([[[ZZ(-7)]]], 2, ZZ) == [[[ZZ(7)]]]

    assert dmp_neg([[[]]], 2, QQ) == [[[]]]
    assert dmp_neg([[[QQ(1,9)]]], 2, QQ) == [[[QQ(-1,9)]]]
    assert dmp_neg([[[QQ(-7,9)]]], 2, QQ) == [[[QQ(7,9)]]]

def test_dup_add():
    assert dup_add([], [], ZZ) == []
    assert dup_add([ZZ(1)], [], ZZ) == [ZZ(1)]
    assert dup_add([], [ZZ(1)], ZZ) == [ZZ(1)]
    assert dup_add([ZZ(1)], [ZZ(1)], ZZ) == [ZZ(2)]
    assert dup_add([ZZ(1)], [ZZ(2)], ZZ) == [ZZ(3)]

    assert dup_add([ZZ(1),ZZ(2)], [ZZ(1)], ZZ) == [ZZ(1),ZZ(3)]
    assert dup_add([ZZ(1)], [ZZ(1),ZZ(2)], ZZ) == [ZZ(1),ZZ(3)]

    assert dup_add([ZZ(1),ZZ(2),ZZ(3)], [ZZ(8),ZZ(9),ZZ(10)], ZZ) == [ZZ(9),ZZ(11),ZZ(13)]

    assert dup_add([], [], QQ) == []
    assert dup_add([QQ(1,2)], [], QQ) == [QQ(1,2)]
    assert dup_add([], [QQ(1,2)], QQ) == [QQ(1,2)]
    assert dup_add([QQ(1,4)], [QQ(1,4)], QQ) == [QQ(1,2)]
    assert dup_add([QQ(1,4)], [QQ(1,2)], QQ) == [QQ(3,4)]

    assert dup_add([QQ(1,2),QQ(2,3)], [QQ(1)], QQ) == [QQ(1,2),QQ(5,3)]
    assert dup_add([QQ(1)], [QQ(1,2),QQ(2,3)], QQ) == [QQ(1,2),QQ(5,3)]

    assert dup_add([QQ(1,7),QQ(2,7),QQ(3,7)], [QQ(8,7),QQ(9,7),QQ(10,7)], QQ) == [QQ(9,7),QQ(11,7),QQ(13,7)]

def test_dmp_add():
    assert dmp_add([ZZ(1),ZZ(2)], [ZZ(1)], 0, ZZ) == \
           dup_add([ZZ(1),ZZ(2)], [ZZ(1)], ZZ)
    assert dmp_add([QQ(1,2),QQ(2,3)], [QQ(1)], 0, QQ) == \
           dup_add([QQ(1,2),QQ(2,3)], [QQ(1)], QQ)

    assert dmp_add([[[]]], [[[]]], 2, ZZ) == [[[]]]
    assert dmp_add([[[ZZ(1)]]], [[[]]], 2, ZZ) == [[[ZZ(1)]]]
    assert dmp_add([[[]]], [[[ZZ(1)]]], 2, ZZ) == [[[ZZ(1)]]]
    assert dmp_add([[[ZZ(2)]]], [[[ZZ(1)]]], 2, ZZ) == [[[ZZ(3)]]]
    assert dmp_add([[[ZZ(1)]]], [[[ZZ(2)]]], 2, ZZ) == [[[ZZ(3)]]]

    assert dmp_add([[[]]], [[[]]], 2, QQ) == [[[]]]
    assert dmp_add([[[QQ(1,2)]]], [[[]]], 2, QQ) == [[[QQ(1,2)]]]
    assert dmp_add([[[]]], [[[QQ(1,2)]]], 2, QQ) == [[[QQ(1,2)]]]
    assert dmp_add([[[QQ(2,7)]]], [[[QQ(1,7)]]], 2, QQ) == [[[QQ(3,7)]]]
    assert dmp_add([[[QQ(1,7)]]], [[[QQ(2,7)]]], 2, QQ) == [[[QQ(3,7)]]]

def test_dup_sub():
    assert dup_sub([], [], ZZ) == []
    assert dup_sub([ZZ(1)], [], ZZ) == [ZZ(1)]
    assert dup_sub([], [ZZ(1)], ZZ) == [ZZ(-1)]
    assert dup_sub([ZZ(1)], [ZZ(1)], ZZ) == []
    assert dup_sub([ZZ(1)], [ZZ(2)], ZZ) == [ZZ(-1)]

    assert dup_sub([ZZ(1),ZZ(2)], [ZZ(1)], ZZ) == [ZZ(1),ZZ(1)]
    assert dup_sub([ZZ(1)], [ZZ(1),ZZ(2)], ZZ) == [ZZ(-1),ZZ(-1)]

    assert dup_sub([ZZ(3),ZZ(2),ZZ(1)], [ZZ(8),ZZ(9),ZZ(10)], ZZ) == [ZZ(-5),ZZ(-7),ZZ(-9)]

    assert dup_sub([], [], QQ) == []
    assert dup_sub([QQ(1,2)], [], QQ) == [QQ(1,2)]
    assert dup_sub([], [QQ(1,2)], QQ) == [QQ(-1,2)]
    assert dup_sub([QQ(1,3)], [QQ(1,3)], QQ) == []
    assert dup_sub([QQ(1,3)], [QQ(2,3)], QQ) == [QQ(-1,3)]

    assert dup_sub([QQ(1,7),QQ(2,7)], [QQ(1)], QQ) == [QQ(1,7),QQ(-5,7)]
    assert dup_sub([QQ(1)], [QQ(1,7),QQ(2,7)], QQ) == [QQ(-1,7),QQ(5,7)]

    assert dup_sub([QQ(3,7),QQ(2,7),QQ(1,7)], [QQ(8,7),QQ(9,7),QQ(10,7)], QQ) == [QQ(-5,7),QQ(-7,7),QQ(-9,7)]

def test_dmp_sub():
    assert dmp_sub([ZZ(1),ZZ(2)], [ZZ(1)], 0, ZZ) == \
           dup_sub([ZZ(1),ZZ(2)], [ZZ(1)], ZZ)
    assert dmp_sub([QQ(1,2),QQ(2,3)], [QQ(1)], 0, QQ) == \
           dup_sub([QQ(1,2),QQ(2,3)], [QQ(1)], QQ)

    assert dmp_sub([[[]]], [[[]]], 2, ZZ) == [[[]]]
    assert dmp_sub([[[ZZ(1)]]], [[[]]], 2, ZZ) == [[[ZZ(1)]]]
    assert dmp_sub([[[]]], [[[ZZ(1)]]], 2, ZZ) == [[[ZZ(-1)]]]
    assert dmp_sub([[[ZZ(2)]]], [[[ZZ(1)]]], 2, ZZ) == [[[ZZ(1)]]]
    assert dmp_sub([[[ZZ(1)]]], [[[ZZ(2)]]], 2, ZZ) == [[[ZZ(-1)]]]

    assert dmp_sub([[[]]], [[[]]], 2, QQ) == [[[]]]
    assert dmp_sub([[[QQ(1,2)]]], [[[]]], 2, QQ) == [[[QQ(1,2)]]]
    assert dmp_sub([[[]]], [[[QQ(1,2)]]], 2, QQ) == [[[QQ(-1,2)]]]
    assert dmp_sub([[[QQ(2,7)]]], [[[QQ(1,7)]]], 2, QQ) == [[[QQ(1,7)]]]
    assert dmp_sub([[[QQ(1,7)]]], [[[QQ(2,7)]]], 2, QQ) == [[[QQ(-1,7)]]]

def test_dup_add_mul():
    assert dup_add_mul([ZZ(1),ZZ(2),ZZ(3)], [ZZ(3),ZZ(2),ZZ(1)],
               [ZZ(1),ZZ(2)], ZZ) == [ZZ(3), ZZ(9), ZZ(7), ZZ(5)]

def test_dup_add_mul():
    assert dmp_add_mul([[ZZ(1),ZZ(2)],[ZZ(3)]], [[ZZ(3)],[ZZ(2),ZZ(1)]],
               [[ZZ(1)],[ZZ(2)]], 1, ZZ) == [[ZZ(3)], [ZZ(3), ZZ(9)], [ZZ(4), ZZ(5)]]

def test_dup_sub_mul():
    assert dup_sub_mul([ZZ(1),ZZ(2),ZZ(3)], [ZZ(3),ZZ(2),ZZ(1)],
               [ZZ(1),ZZ(2)], ZZ) == [ZZ(-3),ZZ(-7),ZZ(-3), ZZ(1)]

def test_dup_sub_mul():
    assert dmp_sub_mul([[ZZ(1),ZZ(2)],[ZZ(3)]], [[ZZ(3)],[ZZ(2),ZZ(1)]],
               [[ZZ(1)],[ZZ(2)]], 1, ZZ) == [[ZZ(-3)], [ZZ(-1), ZZ(-5)], [ZZ(-4), ZZ(1)]]

def test_dup_mul():
    assert dup_mul([], [], ZZ) == []
    assert dup_mul([], [ZZ(1)], ZZ) == []
    assert dup_mul([ZZ(1)], [], ZZ) == []
    assert dup_mul([ZZ(1)], [ZZ(1)], ZZ) == [ZZ(1)]
    assert dup_mul([ZZ(5)], [ZZ(7)], ZZ) == [ZZ(35)]

    assert dup_mul([], [], QQ) == []
    assert dup_mul([], [QQ(1,2)], QQ) == []
    assert dup_mul([QQ(1,2)], [], QQ) == []
    assert dup_mul([QQ(1,2)], [QQ(4,7)], QQ) == [QQ(2,7)]
    assert dup_mul([QQ(5,7)], [QQ(3,7)], QQ) == [QQ(15,49)]

    f = dup_normal([3,0,0,6,1,2], ZZ)
    g = dup_normal([4,0,1,0], ZZ)
    h = dup_normal([12,0,3,24,4,14,1,2,0], ZZ)

    assert dup_mul(f, g, ZZ) == h
    assert dup_mul(g, f, ZZ) == h

    f = dup_normal([2,0,0,1,7], ZZ)
    h = dup_normal([4,0,0,4,28,0,1,14,49], ZZ)

    assert dup_mul(f, f, ZZ) == h

    K = FF(6)

    assert dup_mul([K(2),K(1)], [K(3),K(4)], K) == [K(5),K(4)]

def test_dmp_mul():
    assert dmp_mul([ZZ(5)], [ZZ(7)], 0, ZZ) == \
           dup_mul([ZZ(5)], [ZZ(7)], ZZ)
    assert dmp_mul([QQ(5,7)], [QQ(3,7)], 0, QQ) == \
           dup_mul([QQ(5,7)], [QQ(3,7)], QQ)

    assert dmp_mul([[[]]], [[[]]], 2, ZZ) == [[[]]]
    assert dmp_mul([[[ZZ(1)]]], [[[]]], 2, ZZ) == [[[]]]
    assert dmp_mul([[[]]], [[[ZZ(1)]]], 2, ZZ) == [[[]]]
    assert dmp_mul([[[ZZ(2)]]], [[[ZZ(1)]]], 2, ZZ) == [[[ZZ(2)]]]
    assert dmp_mul([[[ZZ(1)]]], [[[ZZ(2)]]], 2, ZZ) == [[[ZZ(2)]]]

    assert dmp_mul([[[]]], [[[]]], 2, QQ) == [[[]]]
    assert dmp_mul([[[QQ(1,2)]]], [[[]]], 2, QQ) == [[[]]]
    assert dmp_mul([[[]]], [[[QQ(1,2)]]], 2, QQ) == [[[]]]
    assert dmp_mul([[[QQ(2,7)]]], [[[QQ(1,3)]]], 2, QQ) == [[[QQ(2,21)]]]
    assert dmp_mul([[[QQ(1,7)]]], [[[QQ(2,3)]]], 2, QQ) == [[[QQ(2,21)]]]

    K = FF(6)

    assert dmp_mul([[K(2)],[K(1)]], [[K(3)],[K(4)]], 1, K) == [[K(5)],[K(4)]]

def test_dup_sqr():
    assert dup_sqr([], ZZ) == []
    assert dup_sqr([ZZ(2)], ZZ) == [ZZ(4)]
    assert dup_sqr([ZZ(1),ZZ(2)], ZZ) == [ZZ(1),ZZ(4),ZZ(4)]

    assert dup_sqr([], QQ) == []
    assert dup_sqr([QQ(2,3)], QQ) == [QQ(4,9)]
    assert dup_sqr([QQ(1,3),QQ(2,3)], QQ) == [QQ(1,9),QQ(4,9),QQ(4,9)]

    f = dup_normal([2,0,0,1,7], ZZ)

    assert dup_sqr(f, ZZ) == dup_normal([4,0,0,4,28,0,1,14,49], ZZ)

    K = FF(9)

    assert dup_sqr([K(3),K(4)], K) == [K(6),K(7)]

def test_dmp_sqr():
    assert dmp_sqr([ZZ(1),ZZ(2)], 0, ZZ) == \
           dup_sqr([ZZ(1),ZZ(2)], ZZ)

    assert dmp_sqr([[[]]], 2, ZZ) == [[[]]]
    assert dmp_sqr([[[ZZ(2)]]], 2, ZZ) == [[[ZZ(4)]]]

    assert dmp_sqr([[[]]], 2, QQ) == [[[]]]
    assert dmp_sqr([[[QQ(2,3)]]], 2, QQ) == [[[QQ(4,9)]]]

    K = FF(9)

    assert dmp_sqr([[K(3)],[K(4)]], 1, K) == [[K(6)],[K(7)]]

def test_dup_pow():
    assert dup_pow([], 0, ZZ) == [ZZ(1)]
    assert dup_pow([], 0, QQ) == [QQ(1)]

    assert dup_pow([], 1, ZZ) == []
    assert dup_pow([], 7, ZZ) == []

    assert dup_pow([ZZ(1)], 0, ZZ) == [ZZ(1)]
    assert dup_pow([ZZ(1)], 1, ZZ) == [ZZ(1)]
    assert dup_pow([ZZ(1)], 7, ZZ) == [ZZ(1)]

    assert dup_pow([ZZ(3)], 0, ZZ) == [ZZ(1)]
    assert dup_pow([ZZ(3)], 1, ZZ) == [ZZ(3)]
    assert dup_pow([ZZ(3)], 7, ZZ) == [ZZ(2187)]

    assert dup_pow([QQ(1,1)], 0, QQ) == [QQ(1,1)]
    assert dup_pow([QQ(1,1)], 1, QQ) == [QQ(1,1)]
    assert dup_pow([QQ(1,1)], 7, QQ) == [QQ(1,1)]

    assert dup_pow([QQ(3,7)], 0, QQ) == [QQ(1,1)]
    assert dup_pow([QQ(3,7)], 1, QQ) == [QQ(3,7)]
    assert dup_pow([QQ(3,7)], 7, QQ) == [QQ(2187,823543)]

    f = dup_normal([2,0,0,1,7], ZZ)

    assert dup_pow(f, 0, ZZ) == dup_normal([1], ZZ)
    assert dup_pow(f, 1, ZZ) == dup_normal([2,0,0,1,7], ZZ)
    assert dup_pow(f, 2, ZZ) == dup_normal([4,0,0,4,28,0,1,14,49], ZZ)
    assert dup_pow(f, 3, ZZ) == dup_normal([8,0,0,12,84,0,6,84,294,1,21,147,343], ZZ)

def test_dmp_pow():
    assert dmp_pow([[]], 0, 1, ZZ) == [[ZZ(1)]]
    assert dmp_pow([[]], 0, 1, QQ) == [[QQ(1)]]

    assert dmp_pow([[]], 1, 1, ZZ) == [[]]
    assert dmp_pow([[]], 7, 1, ZZ) == [[]]

    assert dmp_pow([[ZZ(1)]], 0, 1, ZZ) == [[ZZ(1)]]
    assert dmp_pow([[ZZ(1)]], 1, 1, ZZ) == [[ZZ(1)]]
    assert dmp_pow([[ZZ(1)]], 7, 1, ZZ) == [[ZZ(1)]]

    assert dmp_pow([[QQ(3,7)]], 0, 1, QQ) == [[QQ(1,1)]]
    assert dmp_pow([[QQ(3,7)]], 1, 1, QQ) == [[QQ(3,7)]]
    assert dmp_pow([[QQ(3,7)]], 7, 1, QQ) == [[QQ(2187,823543)]]

    f = dup_normal([2,0,0,1,7], ZZ)

    assert dmp_pow(f, 2, 0, ZZ) == dup_pow(f, 2, ZZ)

def test_dup_pdiv():
    f = dup_normal([3,1,1,5], ZZ)
    g = dup_normal([5,-3,1], ZZ)

    q = dup_normal([15, 14], ZZ)
    r = dup_normal([52, 111], ZZ)

    assert dup_pdiv(f, g, ZZ) == (q, r)
    assert dup_pexquo(f, g, ZZ) == q
    assert dup_prem(f, g, ZZ) == r

    raises(ExactQuotientFailed, 'dup_pquo(f, g, ZZ)')

    f = dup_normal([3,1,1,5], QQ)
    g = dup_normal([5,-3,1], QQ)

    q = dup_normal([15, 14], QQ)
    r = dup_normal([52, 111], QQ)

    assert dup_pdiv(f, g, QQ) == (q, r)
    assert dup_pexquo(f, g, QQ) == q
    assert dup_prem(f, g, QQ) == r

    raises(ExactQuotientFailed, 'dup_pquo(f, g, QQ)')

def test_dmp_pdiv():
    f = dmp_normal([[1], [], [1,0,0]], 1, ZZ)
    g = dmp_normal([[1], [-1,0]], 1, ZZ)

    q = dmp_normal([[1], [1, 0]], 1, ZZ)
    r = dmp_normal([[2, 0, 0]], 1, ZZ)

    assert dmp_pdiv(f, g, 1, ZZ) == (q, r)
    assert dmp_pexquo(f, g, 1, ZZ) == q
    assert dmp_prem(f, g, 1, ZZ) == r

    raises(ExactQuotientFailed, 'dmp_pquo(f, g, 1, ZZ)')

    f = dmp_normal([[1], [], [1,0,0]], 1, ZZ)
    g = dmp_normal([[2], [-2,0]], 1, ZZ)

    q = dmp_normal([[2], [2, 0]], 1, ZZ)
    r = dmp_normal([[8, 0, 0]], 1, ZZ)

    assert dmp_pdiv(f, g, 1, ZZ) == (q, r)
    assert dmp_pexquo(f, g, 1, ZZ) == q
    assert dmp_prem(f, g, 1, ZZ) == r

    raises(ExactQuotientFailed, 'dmp_pquo(f, g, 1, ZZ)')

def test_dup_rr_div():
    raises(ZeroDivisionError, "dup_rr_div([1,2,3], [], ZZ)")

    f = dup_normal([3,1,1,5], ZZ)
    g = dup_normal([5,-3,1], ZZ)

    q, r = [], f

    assert dup_rr_div(f, g, ZZ) == (q, r)

def test_dmp_rr_div():
    raises(ZeroDivisionError, "dmp_rr_div([[1,2],[3]], [[]], 1, ZZ)")

    f = dmp_normal([[1], [], [1,0,0]], 1, ZZ)
    g = dmp_normal([[1], [-1,0]], 1, ZZ)

    q = dmp_normal([[1], [1, 0]], 1, ZZ)
    r = dmp_normal([[2, 0, 0]], 1, ZZ)

    assert dmp_rr_div(f, g, 1, ZZ) == (q, r)

    f = dmp_normal([[1], [], [1,0,0]], 1, ZZ)
    g = dmp_normal([[-1], [1,0]], 1, ZZ)

    q = dmp_normal([[-1], [-1, 0]], 1, ZZ)
    r = dmp_normal([[2, 0, 0]], 1, ZZ)

    assert dmp_rr_div(f, g, 1, ZZ) == (q, r)

    f = dmp_normal([[1], [], [1,0,0]], 1, ZZ)
    g = dmp_normal([[2], [-2,0]], 1, ZZ)

    q, r = [[]], f

    assert dmp_rr_div(f, g, 1, ZZ) == (q, r)

def test_dup_ff_div():
    raises(ZeroDivisionError, "dup_ff_div([1,2,3], [], QQ)")

    f = dup_normal([3,1,1,5], QQ)
    g = dup_normal([5,-3,1], QQ)

    q = [QQ(3,5), QQ(14,25)]
    r = [QQ(52,25), QQ(111,25)]

    assert dup_ff_div(f, g, QQ) == (q, r)

def test_dmp_ff_div():
    raises(ZeroDivisionError, "dmp_ff_div([[1,2],[3]], [[]], 1, QQ)")

    f = dmp_normal([[1], [], [1,0,0]], 1, QQ)
    g = dmp_normal([[1], [-1,0]], 1, QQ)

    q = [[QQ(1, 1)], [QQ(1, 1), QQ(0, 1)]]
    r = [[QQ(2, 1), QQ(0, 1), QQ(0, 1)]]

    assert dmp_ff_div(f, g, 1, QQ) == (q, r)

    f = dmp_normal([[1], [], [1,0,0]], 1, QQ)
    g = dmp_normal([[-1], [1,0]], 1, QQ)

    q = [[QQ(-1, 1)], [QQ(-1, 1), QQ(0, 1)]]
    r = [[QQ(2, 1), QQ(0, 1), QQ(0, 1)]]

    assert dmp_ff_div(f, g, 1, QQ) == (q, r)

    f = dmp_normal([[1], [], [1,0,0]], 1, QQ)
    g = dmp_normal([[2], [-2,0]], 1, QQ)

    q = [[QQ(1, 2)], [QQ(1, 2), QQ(0, 1)]]
    r = [[QQ(2, 1), QQ(0, 1), QQ(0, 1)]]

    assert dmp_ff_div(f, g, 1, QQ) == (q, r)

def test_dup_div():
    f, g, q, r = [5,4,3,2,1], [1,2,3], [5,-6,0], [20,1]

    assert dup_div(f, g, ZZ) == (q, r)
    assert dup_exquo(f, g, ZZ) == q
    assert dup_rem(f, g, ZZ) == r

    raises(ExactQuotientFailed, 'dup_quo(f, g, ZZ)')

    f, g, q, r = [5,4,3,2,1,0], [1,2,0,0,9], [5,-6], [15,2,-44,54]

    assert dup_div(f, g, ZZ) == (q, r)
    assert dup_exquo(f, g, ZZ) == q
    assert dup_rem(f, g, ZZ) == r

    raises(ExactQuotientFailed, 'dup_quo(f, g, ZZ)')

def test_dmp_div():
    f, g, q, r = [5,4,3,2,1], [1,2,3], [5,-6,0], [20,1]

    assert dmp_div(f, g, 0, ZZ) == (q, r)
    assert dmp_exquo(f, g, 0, ZZ) == q
    assert dmp_rem(f, g, 0, ZZ) == r

    raises(ExactQuotientFailed, 'dmp_quo(f, g, 0, ZZ)')

    f, g, q, r = [[[1]]], [[[2]],[1]], [[[]]], [[[1]]]

    assert dmp_div(f, g, 2, ZZ) == (q, r)
    assert dmp_exquo(f, g, 2, ZZ) == q
    assert dmp_rem(f, g, 2, ZZ) == r

    raises(ExactQuotientFailed, 'dmp_quo(f, g, 2, ZZ)')

def test_dup_max_norm():
    assert dup_max_norm([], ZZ) == 0
    assert dup_max_norm([1], ZZ) == 1

    assert dup_max_norm([1,4,2,3], ZZ) == 4

def test_dmp_max_norm():
    assert dmp_max_norm([[[]]], 2, ZZ) == 0
    assert dmp_max_norm([[[1]]], 2, ZZ) == 1

    assert dmp_max_norm(f_0, 2, ZZ) == 6

def test_dup_l1_norm():
    assert dup_l1_norm([], ZZ) == 0
    assert dup_l1_norm([1], ZZ) == 1
    assert dup_l1_norm([1,4,2,3], ZZ) == 10

def test_dmp_l1_norm():
    assert dmp_l1_norm([[[]]], 2, ZZ) == 0
    assert dmp_l1_norm([[[1]]], 2, ZZ) == 1

    assert dmp_l1_norm(f_0, 2, ZZ) == 31

def test_dup_expand():
    assert dup_expand((), ZZ) == [1]
    assert dup_expand(([1,2,3], [1,2], [7,5,4,3]), ZZ) == \
        dup_mul([1,2,3], dup_mul([1,2], [7,5,4,3], ZZ), ZZ)

def test_dmp_expand():
    assert dmp_expand((), 1, ZZ) == [[1]]
    assert dmp_expand(([[1],[2],[3]], [[1],[2]], [[7],[5],[4],[3]]), 1, ZZ) == \
        dmp_mul([[1],[2],[3]], dmp_mul([[1],[2]], [[7],[5],[4],[3]], 1, ZZ), 1, ZZ)

def test_dup_revert():
    f = [-QQ(1,720),QQ(0),QQ(1,24),QQ(0),-QQ(1,2),QQ(0),QQ(1)]
    g = [QQ(61,720),QQ(0),QQ(5,24),QQ(0), QQ(1,2),QQ(0),QQ(1)]

    assert dup_revert(f, 8, QQ) == g

    raises(NotReversible, "dup_revert([QQ(1), QQ(0)], 3, QQ)")

