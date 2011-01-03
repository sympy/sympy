"""Tests for square--free decomposition algorithms and related tools. """

from sympy.polys.sqfreetools import (
    dup_sqf_p, dmp_sqf_p,
    dup_sqf_norm, dmp_sqf_norm,
    dup_sqf_part, dmp_sqf_part,
    dup_sqf_list, dup_sqf_list_include,
    dmp_sqf_list, dmp_sqf_list_include,
    dup_gff_list, dmp_gff_list)

from sympy.polys.euclidtools import (
    dmp_resultant)

from sympy.polys.densearith import (
    dmp_neg, dmp_sub, dmp_mul, dmp_sqr)

from sympy.polys.densetools import (
    dmp_diff)

from sympy.polys.polyerrors import (
    DomainError)

from sympy.polys.specialpolys import (
    f_0, f_1, f_2, f_3, f_4, f_5, f_6)

from sympy.polys.domains import FF, ZZ, QQ

from sympy.utilities.pytest import raises

def test_dup_sqf():
    assert dup_sqf_part([], ZZ) == []
    assert dup_sqf_p([], ZZ) == True

    assert dup_sqf_part([7], ZZ) == [1]
    assert dup_sqf_p([7], ZZ) == True

    assert dup_sqf_part([2,2], ZZ) == [1,1]
    assert dup_sqf_p([2,2], ZZ) == True

    assert dup_sqf_part([1,0,1,1], ZZ) == [1,0,1,1]
    assert dup_sqf_p([1,0,1,1], ZZ) == True

    assert dup_sqf_part([-1,0,1,1], ZZ) == [1,0,-1,-1]
    assert dup_sqf_p([-1,0,1,1], ZZ) == True

    assert dup_sqf_part([2,3,0,0], ZZ) == [2,3,0]
    assert dup_sqf_p([2,3,0,0], ZZ) == False

    assert dup_sqf_part([-2,3,0,0], ZZ) == [2,-3,0]
    assert dup_sqf_p([-2,3,0,0], ZZ) == False

    assert dup_sqf_list([], ZZ) == (0, [])
    assert dup_sqf_list([1], ZZ) == (1, [])

    assert dup_sqf_list([1,0], ZZ) == (1, [([1,0], 1)])
    assert dup_sqf_list([2,0,0], ZZ) == (2, [([1,0], 2)])
    assert dup_sqf_list([3,0,0,0], ZZ) == (3, [([1,0], 3)])

    assert dup_sqf_list([ZZ(2),ZZ(4),ZZ(2)], ZZ) == \
        (ZZ(2), [([ZZ(1),ZZ(1)], 2)])
    assert dup_sqf_list([QQ(2),QQ(4),QQ(2)], QQ) == \
        (QQ(2), [([QQ(1),QQ(1)], 2)])

    assert dup_sqf_list([-1,1,0,0,1,-1], ZZ) == \
        (-1, [([1,1,1,1], 1), ([1,-1], 2)])
    assert dup_sqf_list([1,0,6,0,12,0,8,0,0], ZZ) == \
        (1, [([1,0], 2), ([1,0,2], 3)])

    K = FF(2)
    f = map(K, [1,0,1])

    assert dup_sqf_list(f, K) == \
        (K(1), [([K(1),K(1)], 2)])

    K = FF(3)
    f = map(K, [1,0,0,2,0,0,2,0,0,1,0])

    assert dup_sqf_list(f, K) == \
        (K(1), [([K(1), K(0)], 1),
                ([K(1), K(1)], 3),
                ([K(1), K(2)], 6)])

    f = [1,0,0,1]
    g = map(K, f)

    assert dup_sqf_part(f, ZZ) == f
    assert dup_sqf_part(g, K) == [K(1), K(1)]

    assert dup_sqf_p(f, ZZ) == True
    assert dup_sqf_p(g, K) == False

    A = [[1],[],[-3],[],[6]]
    D = [[1],[],[-5],[],[5],[],[4]]

    f, g = D, dmp_sub(A, dmp_mul(dmp_diff(D, 1, 1, ZZ), [[1,0]], 1, ZZ), 1, ZZ)

    res = dmp_resultant(f, g, 1, ZZ)

    assert dup_sqf_list(res, ZZ) == (45796, [([4,0,1], 3)])

def test_dmp_sqf():
    assert dmp_sqf_part([[]], 1, ZZ) == [[]]
    assert dmp_sqf_p([[]], 1, ZZ) == True

    assert dmp_sqf_part([[7]], 1, ZZ) == [[1]]
    assert dmp_sqf_p([[7]], 1, ZZ) == True

    assert dmp_sqf_p(f_0, 2, ZZ) == True
    assert dmp_sqf_p(dmp_sqr(f_0, 2, ZZ), 2, ZZ) == False
    assert dmp_sqf_p(f_1, 2, ZZ) == True
    assert dmp_sqf_p(dmp_sqr(f_1, 2, ZZ), 2, ZZ) == False
    assert dmp_sqf_p(f_2, 2, ZZ) == True
    assert dmp_sqf_p(dmp_sqr(f_2, 2, ZZ), 2, ZZ) == False
    assert dmp_sqf_p(f_3, 2, ZZ) == True
    assert dmp_sqf_p(dmp_sqr(f_3, 2, ZZ), 2, ZZ) == False
    assert dmp_sqf_p(f_5, 2, ZZ) == False
    assert dmp_sqf_p(dmp_sqr(f_5, 2, ZZ), 2, ZZ) == False

    assert dmp_sqf_p(f_4, 2, ZZ) == True
    assert dmp_sqf_part(f_4, 2, ZZ) == dmp_neg(f_4, 2, ZZ)
    assert dmp_sqf_p(f_6, 3, ZZ) == True
    assert dmp_sqf_part(f_6, 3, ZZ) == f_6

    assert dmp_sqf_part(f_5, 2, ZZ) == [[[1]], [[1], [-1, 0]]]

    assert dup_sqf_list([], ZZ) == (ZZ(0), [])
    assert dup_sqf_list_include([], ZZ) == [([], 1)]

    assert dmp_sqf_list([[ZZ(3)]], 1, ZZ) == (ZZ(3), [])
    assert dmp_sqf_list_include([[ZZ(3)]], 1, ZZ) == [([[ZZ(3)]], 1)]

    f = [-1,1,0,0,1,-1]

    assert dmp_sqf_list(f, 0, ZZ) == \
        (-1, [([1,1,1,1], 1), ([1,-1], 2)])
    assert dmp_sqf_list_include(f, 0, ZZ) == \
        [([-1,-1,-1,-1], 1), ([1,-1], 2)]

    f = [[-1],[1],[],[],[1],[-1]]

    assert dmp_sqf_list(f, 1, ZZ) == \
        (-1, [([[1],[1],[1],[1]], 1), ([[1],[-1]], 2)])
    assert dmp_sqf_list_include(f, 1, ZZ) == \
        [([[-1],[-1],[-1],[-1]], 1), ([[1],[-1]], 2)]

    K = FF(2)

    raises(DomainError, "dmp_sqf_list([[K(1), K(0), K(1)]], 1, K)")

def test_dup_gff_list():
    f = [1, 2, -1, -2, 0, 0]

    assert dup_gff_list(f, ZZ) == [([1, 0], 1), ([1, 2], 4)]

    g = [1, -20, 166, -744, 1965, -3132, 2948, -1504, 320, 0]

    assert dup_gff_list(g, ZZ) == [([1, -5, 4], 1), ([1, -5, 4], 2), ([1, 0], 3)]

    raises(ValueError, "dup_gff_list([], ZZ)")
