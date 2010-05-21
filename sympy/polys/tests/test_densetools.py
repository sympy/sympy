"""Tests for dense recursive polynomials' tools. """

from sympy.polys.densebasic import (
    dup_LC, dmp_LC, dup_normal, dmp_normal,
    dup_from_raw_dict, dmp_from_dict,
    dmp_convert, dmp_swap, dmp_one_p,
)

from sympy.polys.densearith import (
    dup_add, dup_mul, dup_exquo,
    dmp_neg, dmp_sub, dmp_mul_ground, dmp_mul, dmp_sqr,
)

from sympy.polys.densetools import (
    dup_clear_denoms, dmp_clear_denoms,
    dup_integrate, dmp_integrate, dmp_integrate_in,
    dup_diff, dmp_diff, dmp_diff_in,
    dup_eval, dmp_eval, dmp_eval_in,
    dmp_eval_tail, dmp_diff_eval_in,
    dup_gcdex, dup_half_gcdex, dup_invert,
    dup_subresultants, dmp_subresultants,
    dup_prs_resultant, dmp_prs_resultant,
    dmp_zz_collins_resultant,
    dmp_qq_collins_resultant,
    dup_resultant, dmp_resultant,
    dup_discriminant, dmp_discriminant,
    dup_rr_prs_gcd, dmp_rr_prs_gcd,
    dup_ff_prs_gcd, dmp_ff_prs_gcd,
    dup_zz_heu_gcd, dmp_zz_heu_gcd,
    dup_qq_heu_gcd, dmp_qq_heu_gcd,
    dup_inner_gcd, dmp_inner_gcd,
    dup_gcd, dmp_gcd,
    dup_rr_lcm, dup_ff_lcm, dup_lcm,
    dmp_rr_lcm, dmp_ff_lcm, dmp_lcm,
    dup_trunc, dmp_trunc, dmp_ground_trunc,
    dup_monic, dmp_ground_monic,
    dup_rr_content, dup_ff_content, dup_content,
    dmp_content,
    dmp_rr_ground_content, dmp_ff_ground_content, dmp_ground_content,
    dup_rr_primitive, dup_ff_primitive, dup_primitive,
    dmp_primitive,
    dmp_rr_ground_primitive, dmp_ff_ground_primitive, dmp_ground_primitive,
    dup_sqf_p, dmp_sqf_p,
    dup_sqf_part, dmp_sqf_part,
    dup_sqf_list, dup_sqf_list_include,
    dmp_sqf_list, dmp_sqf_list_include,
    dup_extract, dmp_ground_extract,
    dup_real_imag,
    dup_mirror, dup_scale, dup_taylor,
    dup_transform,
    dup_compose, dmp_compose,
    dup_decompose,
    dup_sturm,
    dmp_lift,
    dup_sign_variations,
)

from sympy.polys.polyclasses import DMP, ANP

from sympy.polys.polyerrors import (
    ExactQuotientFailed,
    HeuristicGCDFailed,
    RefinementFailed,
    DomainError,
)

from sympy.polys.specialpolys import (
    f_0, f_1, f_2, f_3, f_4, f_5, f_6,
    dmp_fateman_poly_F_1,
    dmp_fateman_poly_F_2,
    dmp_fateman_poly_F_3,
)

from sympy.polys.algebratools import ZZ, QQ, EX

from sympy import raises, I

def test_dup_clear_denoms():
    assert dup_clear_denoms([], QQ, ZZ) == (ZZ(1), [])

    assert dup_clear_denoms([QQ(1)], QQ, ZZ) == (ZZ(1), [QQ(1)])
    assert dup_clear_denoms([QQ(7)], QQ, ZZ) == (ZZ(1), [QQ(7)])

    assert dup_clear_denoms([QQ(7,3)], QQ) == (ZZ(3), [QQ(7)])
    assert dup_clear_denoms([QQ(7,3)], QQ, ZZ) == (ZZ(3), [QQ(7)])

    assert dup_clear_denoms([QQ(3),QQ(1),QQ(0)], QQ, ZZ) == (ZZ(1), [QQ(3),QQ(1),QQ(0)])
    assert dup_clear_denoms([QQ(1),QQ(1,2),QQ(0)], QQ, ZZ) == (ZZ(2), [QQ(2),QQ(1),QQ(0)])

    assert dup_clear_denoms([QQ(3),QQ(1),QQ(0)], QQ, ZZ, convert=True) == (ZZ(1), [ZZ(3),ZZ(1),ZZ(0)])
    assert dup_clear_denoms([QQ(1),QQ(1,2),QQ(0)], QQ, ZZ, convert=True) == (ZZ(2), [ZZ(2),ZZ(1),ZZ(0)])

    raises(DomainError, "dup_clear_denoms([EX(7)], EX)")

def test_dmp_clear_denoms():
    assert dmp_clear_denoms([[]], 1, QQ, ZZ) == (ZZ(1), [[]])

    assert dmp_clear_denoms([[QQ(1)]], 1, QQ, ZZ) == (ZZ(1), [[QQ(1)]])
    assert dmp_clear_denoms([[QQ(7)]], 1, QQ, ZZ) == (ZZ(1), [[QQ(7)]])

    assert dmp_clear_denoms([[QQ(7,3)]], 1, QQ) == (ZZ(3), [[QQ(7)]])
    assert dmp_clear_denoms([[QQ(7,3)]], 1, QQ, ZZ) == (ZZ(3), [[QQ(7)]])

    assert dmp_clear_denoms([[QQ(3)],[QQ(1)],[]], 1, QQ, ZZ) == (ZZ(1), [[QQ(3)],[QQ(1)],[]])
    assert dmp_clear_denoms([[QQ(1)],[QQ(1,2)],[]], 1, QQ, ZZ) == (ZZ(2), [[QQ(2)],[QQ(1)],[]])

    assert dmp_clear_denoms([[QQ(3)],[QQ(1)],[]], 1, QQ, ZZ, convert=True) == (ZZ(1), [[QQ(3)],[QQ(1)],[]])
    assert dmp_clear_denoms([[QQ(1)],[QQ(1,2)],[]], 1, QQ, ZZ, convert=True) == (ZZ(2), [[QQ(2)],[QQ(1)],[]])

    raises(DomainError, "dmp_clear_denoms([[EX(7)]], 1, EX)")

def test_dup_integrate():
    assert dup_integrate([], 1, QQ) == []
    assert dup_integrate([], 2, QQ) == []

    assert dup_integrate([QQ(1)], 1, QQ) == [QQ(1),QQ(0)]
    assert dup_integrate([QQ(1)], 2, QQ) == [QQ(1,2),QQ(0),QQ(0)]

    assert dup_integrate([QQ(1),QQ(2),QQ(3)], 0, QQ) == \
        [QQ(1),QQ(2),QQ(3)]
    assert dup_integrate([QQ(1),QQ(2),QQ(3)], 1, QQ) == \
        [QQ(1,3),QQ(1),QQ(3),QQ(0)]
    assert dup_integrate([QQ(1),QQ(2),QQ(3)], 2, QQ) == \
        [QQ(1,12),QQ(1,3),QQ(3,2),QQ(0),QQ(0)]
    assert dup_integrate([QQ(1),QQ(2),QQ(3)], 3, QQ) == \
        [QQ(1,60),QQ(1,12),QQ(1,2),QQ(0),QQ(0),QQ(0)]

    assert dup_integrate(dup_from_raw_dict({29: QQ(17)}, QQ), 3, QQ) == \
        dup_from_raw_dict({32: QQ(17,29760)}, QQ)

    assert dup_integrate(dup_from_raw_dict({29: QQ(17), 5: QQ(1,2)}, QQ), 3, QQ) == \
        dup_from_raw_dict({32: QQ(17,29760), 8: QQ(1, 672)}, QQ)

def test_dmp_integrate():
    assert dmp_integrate([[[]]], 1, 2, QQ) == [[[]]]
    assert dmp_integrate([[[]]], 2, 2, QQ) == [[[]]]

    assert dmp_integrate([[[QQ(1)]]], 1, 2, QQ) == [[[QQ(1)]],[[]]]
    assert dmp_integrate([[[QQ(1)]]], 2, 2, QQ) == [[[QQ(1,2)]],[[]],[[]]]

    assert dmp_integrate([[QQ(1)],[QQ(2)],[QQ(3)]], 0, 1, QQ) == \
        [[QQ(1)],[QQ(2)],[QQ(3)]]
    assert dmp_integrate([[QQ(1)],[QQ(2)],[QQ(3)]], 1, 1, QQ) == \
        [[QQ(1,3)],[QQ(1)],[QQ(3)],[]]
    assert dmp_integrate([[QQ(1)],[QQ(2)],[QQ(3)]], 2, 1, QQ) == \
        [[QQ(1,12)],[QQ(1,3)],[QQ(3,2)],[],[]]
    assert dmp_integrate([[QQ(1)],[QQ(2)],[QQ(3)]], 3, 1, QQ) == \
        [[QQ(1,60)],[QQ(1,12)],[QQ(1,2)],[],[],[]]

def test_dmp_integrate_in():
    f = dmp_convert(f_6, 3, ZZ, QQ)

    assert dmp_integrate_in(f, 2, 1, 3, QQ) == \
        dmp_swap(dmp_integrate(dmp_swap(f, 0, 1, 3, QQ), 2, 3, QQ), 0, 1, 3, QQ)
    assert dmp_integrate_in(f, 3, 1, 3, QQ) == \
        dmp_swap(dmp_integrate(dmp_swap(f, 0, 1, 3, QQ), 3, 3, QQ), 0, 1, 3, QQ)
    assert dmp_integrate_in(f, 2, 2, 3, QQ) == \
        dmp_swap(dmp_integrate(dmp_swap(f, 0, 2, 3, QQ), 2, 3, QQ), 0, 2, 3, QQ)
    assert dmp_integrate_in(f, 3, 2, 3, QQ) == \
        dmp_swap(dmp_integrate(dmp_swap(f, 0, 2, 3, QQ), 3, 3, QQ), 0, 2, 3, QQ)

def test_dup_diff():
    assert dup_diff([], 1, ZZ) == []
    assert dup_diff([7], 1, ZZ) == []
    assert dup_diff([2,7], 1, ZZ) == [2]
    assert dup_diff([1,2,1], 1, ZZ) == [2,2]
    assert dup_diff([1,2,3,4], 1, ZZ) == [3,4,3]
    assert dup_diff([1,-1,0,0,2], 1, ZZ) == [4,-3,0,0]

    f = dup_normal([17,34,56,-345,23,76,0,0,12,3,7], ZZ)

    assert dup_diff(f, 0, ZZ) ==                            f
    assert dup_diff(f, 1, ZZ) ==                   dup_diff(f, 1, ZZ)
    assert dup_diff(f, 2, ZZ) ==          dup_diff(dup_diff(f, 1, ZZ), 1, ZZ)
    assert dup_diff(f, 3, ZZ) == dup_diff(dup_diff(dup_diff(f, 1, ZZ), 1, ZZ), 1, ZZ)

def test_dmp_diff():
    assert dmp_diff([], 1, 0, ZZ) == []
    assert dmp_diff([[]], 1, 1, ZZ) == [[]]
    assert dmp_diff([[[]]], 1, 2, ZZ) == [[[]]]

    assert dmp_diff([[[1], [2]]], 1, 2, ZZ) == [[[]]]

    assert dmp_diff([[[1]], [[]]], 1, 2, ZZ) == [[[1]]]
    assert dmp_diff([[[3]], [[1]], [[]]], 1, 2, ZZ) == [[[6]], [[1]]]

    assert dmp_diff([1,-1,0,0,2], 1, 0, ZZ) == \
           dup_diff([1,-1,0,0,2], 1, ZZ)

    assert dmp_diff(f_6, 0, 3, ZZ) ==                            f_6
    assert dmp_diff(f_6, 1, 3, ZZ) ==                   dmp_diff(f_6, 1, 3, ZZ)
    assert dmp_diff(f_6, 2, 3, ZZ) ==          dmp_diff(dmp_diff(f_6, 1, 3, ZZ), 1, 3, ZZ)
    assert dmp_diff(f_6, 3, 3, ZZ) == dmp_diff(dmp_diff(dmp_diff(f_6, 1, 3, ZZ), 1, 3, ZZ), 1, 3, ZZ)

def test_dmp_diff_in():
    assert dmp_diff_in(f_6, 2, 1, 3, ZZ) == \
        dmp_swap(dmp_diff(dmp_swap(f_6, 0, 1, 3, ZZ), 2, 3, ZZ), 0, 1, 3, ZZ)
    assert dmp_diff_in(f_6, 3, 1, 3, ZZ) == \
        dmp_swap(dmp_diff(dmp_swap(f_6, 0, 1, 3, ZZ), 3, 3, ZZ), 0, 1, 3, ZZ)
    assert dmp_diff_in(f_6, 2, 2, 3, ZZ) == \
        dmp_swap(dmp_diff(dmp_swap(f_6, 0, 2, 3, ZZ), 2, 3, ZZ), 0, 2, 3, ZZ)
    assert dmp_diff_in(f_6, 3, 2, 3, ZZ) == \
        dmp_swap(dmp_diff(dmp_swap(f_6, 0, 2, 3, ZZ), 3, 3, ZZ), 0, 2, 3, ZZ)

def test_dup_eval():
    assert dup_eval([], 7, ZZ) == 0
    assert dup_eval([1,2], 0, ZZ) == 2
    assert dup_eval([1,2,3], 7, ZZ) == 66

def test_dmp_eval():
    assert dmp_eval([], 3, 0, ZZ) == 0

    assert dmp_eval([[]], 3, 1, ZZ) == []
    assert dmp_eval([[[]]], 3, 2, ZZ) == [[]]

    assert dmp_eval([[1,2]], 0, 1, ZZ) == [1,2]

    assert dmp_eval([[[1]]], 3, 2, ZZ) == [[1]]
    assert dmp_eval([[[1, 2]]], 3, 2, ZZ) == [[1, 2]]

    assert dmp_eval([[3, 2], [1, 2]], 3, 1, ZZ) == [10, 8]
    assert dmp_eval([[[3, 2]], [[1, 2]]], 3, 2, ZZ) == [[10, 8]]

def test_dmp_eval_in():
    assert dmp_eval_in(f_6,-2, 1, 3, ZZ) == dmp_eval(dmp_swap(f_6, 0, 1, 3, ZZ),-2, 3, ZZ)
    assert dmp_eval_in(f_6, 7, 1, 3, ZZ) == dmp_eval(dmp_swap(f_6, 0, 1, 3, ZZ), 7, 3, ZZ)
    assert dmp_eval_in(f_6,-2, 2, 3, ZZ) == dmp_swap(dmp_eval(dmp_swap(f_6, 0, 2, 3, ZZ),-2, 3, ZZ), 0, 1, 2, ZZ)
    assert dmp_eval_in(f_6, 7, 2, 3, ZZ) == dmp_swap(dmp_eval(dmp_swap(f_6, 0, 2, 3, ZZ), 7, 3, ZZ), 0, 1, 2, ZZ)

    f = [[[45L]], [[]], [[]], [[-9L], [-1L], [], [3L, 0L, 10L, 0L]]]

    assert dmp_eval_in(f, -2, 2, 2, ZZ) == \
        [[45], [], [], [-9, -1, 0, -44]]

def test_dmp_eval_tail():
    assert dmp_eval_tail([[]], [1], 1, ZZ) == []
    assert dmp_eval_tail([[[]]], [1], 2, ZZ) == [[]]
    assert dmp_eval_tail([[[]]], [1, 2], 2, ZZ) == []

    assert dmp_eval_tail(f_0, [], 2, ZZ) == f_0

    assert dmp_eval_tail(f_0, [1,-17,8], 2, ZZ) == 84496
    assert dmp_eval_tail(f_0, [-17, 8], 2, ZZ) == [-1409, 3, 85902]
    assert dmp_eval_tail(f_0, [8], 2, ZZ) == [[83, 2], [3], [302, 81, 1]]

    assert dmp_eval_tail(f_1, [-17, 8], 2, ZZ) == [-136, 15699, 9166, -27144]

    assert dmp_eval_tail(f_2, [-12, 3], 2, ZZ) == [-1377, 0, -702, -1224, 0, -624]
    assert dmp_eval_tail(f_3, [-12, 3], 2, ZZ) == [144, 82, -5181, -28872, -14868, -540]

    assert dmp_eval_tail(f_4, [25, -1], 2, ZZ) == [152587890625, 9765625, -59605407714843750,
        -3839159765625, -1562475, 9536712644531250, 610349546750, -4, 24414375000, 1562520]
    assert dmp_eval_tail(f_5, [25, -1], 2, ZZ) == [-1, -78, -2028, -17576]

    assert dmp_eval_tail(f_6, [0, 2, 4], 3, ZZ) == [5040, 0, 0, 4480]

def test_dmp_diff_eval_in():
    assert dmp_diff_eval_in(f_6, 2, 7, 1, 3, ZZ) == \
        dmp_eval(dmp_diff(dmp_swap(f_6, 0, 1, 3, ZZ), 2, 3, ZZ), 7, 3, ZZ)

def test_dup_gcdex():
    f = dup_normal([1,-2,-6,12,15], QQ)
    g = dup_normal([1,1,-4,-4], QQ)

    s = [QQ(-1,5),QQ(3,5)]
    t = [QQ(1,5),QQ(-6,5),QQ(2)]
    h = [QQ(1),QQ(1)]

    assert dup_half_gcdex(f, g, QQ) == (s, h)
    assert dup_gcdex(f, g, QQ) == (s, t, h)

    f = dup_normal([1,4,0,-1,1], QQ)
    g = dup_normal([1,0,-1,1], QQ)

    s, t, h = dup_gcdex(f, g, QQ)
    S, T, H = dup_gcdex(g, f, QQ)

    assert dup_add(dup_mul(s, f, QQ),
                   dup_mul(t, g, QQ), QQ) == h
    assert dup_add(dup_mul(S, g, QQ),
                   dup_mul(T, f, QQ), QQ) == H

    f = dup_normal([2,0], QQ)
    g = dup_normal([1,0,-16], QQ)

    s = [QQ(1,32),QQ(0)]
    t = [QQ(-1,16)]
    h = [QQ(1)]

    assert dup_half_gcdex(f, g, QQ) == (s, h)
    assert dup_gcdex(f, g, QQ) == (s, t, h)

def test_dup_invert():
    assert dup_invert([QQ(2),QQ(0)], [QQ(1),QQ(0),QQ(-16)], QQ) == [QQ(1,32),QQ(0)]

def test_dup_subresultants():
    assert dup_resultant([], [], ZZ) == ZZ(0)

    assert dup_resultant([ZZ(1)], [], ZZ) == ZZ(0)
    assert dup_resultant([], [ZZ(1)], ZZ) == ZZ(0)

    f = dup_normal([1,0,1,0,-3,-3,8,2,-5], ZZ)
    g = dup_normal([3,0,5,0,-4,-9,21], ZZ)

    a = dup_normal([15,0,-3,0,9], ZZ)
    b = dup_normal([65,125,-245], ZZ)
    c = dup_normal([9326,-12300], ZZ)
    d = dup_normal([260708], ZZ)

    assert dup_subresultants(f, g, ZZ) == [f, g, a, b, c, d]
    assert dup_resultant(f, g, ZZ) == dup_LC(d, ZZ)

    f = dup_normal([1,-2,1], ZZ)
    g = dup_normal([1,0,-1], ZZ)

    a = dup_normal([2,-2], ZZ)

    assert dup_subresultants(f, g, ZZ) == [f, g, a]
    assert dup_resultant(f, g, ZZ) == 0

    f = dup_normal([1,0, 1], ZZ)
    g = dup_normal([1,0,-1], ZZ)

    a = dup_normal([-2], ZZ)

    assert dup_subresultants(f, g, ZZ) ==  [f, g, a]
    assert dup_resultant(f, g, ZZ) == 4

    f = dup_normal([1,0,-1], ZZ)
    g = dup_normal([1,-1,0,2], ZZ)

    assert dup_resultant(f, g, ZZ) == 0

    f = dup_normal([3,0,-1,0], ZZ)
    g = dup_normal([5,0,1], ZZ)

    assert dup_resultant(f, g, ZZ) == 64

    f = dup_normal([1,-2,7], ZZ)
    g = dup_normal([1,0,-1,5], ZZ)

    assert dup_resultant(f, g, ZZ) == 265

    f = dup_normal([1,-6,11,-6], ZZ)
    g = dup_normal([1,-15,74,-120], ZZ)

    assert dup_resultant(f, g, ZZ) == -8640

    f = dup_normal([1,-6,11,-6], ZZ)
    g = dup_normal([1,-10,29,-20], ZZ)

    assert dup_resultant(f, g, ZZ) == 0

    f = dup_normal([1,0,0,-1], ZZ)
    g = dup_normal([1,2,2,-1], ZZ)

    assert dup_resultant(f, g, ZZ) == 16

    f = dup_normal([1,0,0,0,0,0,0,0,-2], ZZ)
    g = dup_normal([1,-1], ZZ)

    assert dup_resultant(f, g, ZZ) == -1

def test_dmp_subresultants():
    assert dmp_resultant([[]], [[]], 1, ZZ) == []
    assert dmp_prs_resultant([[]], [[]], 1, ZZ)[0] == []
    assert dmp_zz_collins_resultant([[]], [[]], 1, ZZ) == []
    assert dmp_qq_collins_resultant([[]], [[]], 1, ZZ) == []

    assert dmp_resultant([[ZZ(1)]], [[]], 1, ZZ) == []
    assert dmp_resultant([[ZZ(1)]], [[]], 1, ZZ) == []
    assert dmp_resultant([[ZZ(1)]], [[]], 1, ZZ) == []

    assert dmp_resultant([[]], [[ZZ(1)]], 1, ZZ) == []
    assert dmp_prs_resultant([[]], [[ZZ(1)]], 1, ZZ)[0] == []
    assert dmp_zz_collins_resultant([[]], [[ZZ(1)]], 1, ZZ) == []
    assert dmp_qq_collins_resultant([[]], [[ZZ(1)]], 1, ZZ) == []

    f = dmp_normal([[3,0],[],[-1,0,0,-4]], 1, ZZ)
    g = dmp_normal([[1],[1,0,0,0],[-9]], 1, ZZ)

    a = dmp_normal([[3,0,0,0,0],[1,0,-27,4]], 1, ZZ)
    b = dmp_normal([[-3,0,0,-12,1,0,-54,8,729,-216,16]], 1, ZZ)

    r = dmp_LC(b, ZZ)

    assert dmp_subresultants(f, g, 1, ZZ) == [f, g, a, b]

    assert dmp_resultant(f, g, 1, ZZ) == r
    assert dmp_prs_resultant(f, g, 1, ZZ)[0] == r
    assert dmp_zz_collins_resultant(f, g, 1, ZZ) == r
    assert dmp_qq_collins_resultant(f, g, 1, ZZ) == r

    f = dmp_normal([[-1],[],[],[5]], 1, ZZ)
    g = dmp_normal([[3,1],[],[]], 1, ZZ)

    a = dmp_normal([[45,30,5]], 1, ZZ)
    b = dmp_normal([[675,675,225,25]], 1, ZZ)

    r = dmp_LC(b, ZZ)

    assert dmp_subresultants(f, g, 1, ZZ) == [f, g, a]
    assert dmp_resultant(f, g, 1, ZZ) == r
    assert dmp_prs_resultant(f, g, 1, ZZ)[0] == r
    assert dmp_zz_collins_resultant(f, g, 1, ZZ) == r
    assert dmp_qq_collins_resultant(f, g, 1, ZZ) == r

    f = [[[[[6]]]], [[[[-3]]], [[[-2]], [[]]]], [[[[1]], [[]]], [[[]]]]]
    g = [[[[[1]]]], [[[[-1], [-1, 0]]]], [[[[1, 0], []]]]]

    r = [[[[1]], [[-3], [-3, 0]], [[9, 0], []]], [[[-2], [-2, 0]], [[6],
         [12, 0], [6, 0, 0]], [[-18, 0], [-18, 0, 0], []]], [[[4, 0],
         []], [[-12, 0], [-12, 0, 0], []], [[36, 0, 0], [], []]]]

    assert dmp_zz_collins_resultant(f, g, 4, ZZ) == r

    f = [[[[[QQ(1,1)]]]], [[[[QQ(-1,2)]]], [[[QQ(-1,3)]], [[]]]], [[[[QQ(1,6)]], [[]]], [[[]]]]]
    g = [[[[[QQ(1,1)]]]], [[[[QQ(-1,1)], [QQ(-1,1), QQ(0, 1)]]]], [[[[QQ(1,1), QQ(0,1)], []]]]]

    r = [[[[QQ(1,36)]], [[QQ(-1,12)], [QQ(-1,12), QQ(0,1)]], [[QQ(1,4), QQ(0,1)], []]],
         [[[QQ(-1,18)], [QQ(-1,18), QQ(0,1)]], [[QQ(1,6)], [QQ(1,3), QQ(0,1)], [QQ(1,6),
            QQ(0,1), QQ(0,1)]], [[QQ(-1,2), QQ(0,1)], [QQ(-1,2), QQ(0,1), QQ(0,1)], []]],
         [[[QQ(1,9), QQ(0,1)], []], [[QQ(-1,3), QQ(0,1)], [QQ(-1,3), QQ(0,1), QQ(0,1)], []],
          [[QQ(1,1), QQ(0,1), QQ(0,1)], [], []]]]

    assert dmp_qq_collins_resultant(f, g, 4, QQ) == r

def test_dup_discriminant():
    assert dup_discriminant([], ZZ) == 0
    assert dup_discriminant([1,0], ZZ) == 1

    assert dup_discriminant([1,3,9,-13], ZZ) == -11664
    assert dup_discriminant([5,0,1,0,0,2], ZZ) == 31252160
    assert dup_discriminant([1,2,6,-22,13], ZZ) == 0
    assert dup_discriminant([12,0,0,15,30,1,0,1], ZZ) == -220289699947514112

def test_dmp_discriminant():
    assert dmp_discriminant([], 0, ZZ) == 0
    assert dmp_discriminant([[]], 1, ZZ) == []

    assert dmp_discriminant([[1,0]], 1, ZZ) == []

    assert dmp_discriminant([1,3,9,-13], 0, ZZ) == -11664
    assert dmp_discriminant([5,0,1,0,0,2], 0, ZZ) == 31252160
    assert dmp_discriminant([1,2,6,-22,13], 0, ZZ) == 0
    assert dmp_discriminant([12,0,0,15,30,1,0,1], 0, ZZ) == -220289699947514112

    assert dmp_discriminant([[1,0],[],[2,0]], 1, ZZ) == [-8,0,0]
    assert dmp_discriminant([[1,0,2],[]], 1, ZZ) == [1]

    assert dmp_discriminant([[[1],[]],[[1,0]]], 2, ZZ) == [[1]]

    assert dmp_discriminant([[[[1]],[[]]],[[[1],[]]],[[[1,0]]]], 3, ZZ) == \
        [[[-4, 0]], [[1], [], []]]
    assert dmp_discriminant([[[[[1]]],[[[]]]],[[[[1]],[[]]]],[[[[1],[]]]],[[[[1,0]]]]], 4, ZZ) == \
        [[[[-27,0,0]]],[[[18,0],[]],[[-4],[],[],[]]],[[[-4,0]],[[1],[],[]],[[]],[[]]]]

def test_dup_gcd():
    assert dup_zz_heu_gcd([], [], ZZ) == ([], [], [])
    assert dup_rr_prs_gcd([], [], ZZ) == ([], [], [])

    assert dup_zz_heu_gcd([2], [], ZZ) == ([2], [1], [])
    assert dup_rr_prs_gcd([2], [], ZZ) == ([2], [1], [])

    assert dup_zz_heu_gcd([-2], [], ZZ) == ([2], [-1], [])
    assert dup_rr_prs_gcd([-2], [], ZZ) == ([2], [-1], [])

    assert dup_zz_heu_gcd([], [-2], ZZ) == ([2], [], [-1])
    assert dup_rr_prs_gcd([], [-2], ZZ) == ([2], [], [-1])

    assert dup_zz_heu_gcd([], [2,4], ZZ) == ([2,4], [], [1])
    assert dup_rr_prs_gcd([], [2,4], ZZ) == ([2,4], [], [1])

    assert dup_zz_heu_gcd([2,4], [], ZZ) == ([2,4], [1], [])
    assert dup_rr_prs_gcd([2,4], [], ZZ) == ([2,4], [1], [])

    assert dup_zz_heu_gcd([2], [2], ZZ) == ([2], [1], [1])
    assert dup_rr_prs_gcd([2], [2], ZZ) == ([2], [1], [1])

    assert dup_zz_heu_gcd([-2], [2], ZZ) == ([2], [-1], [1])
    assert dup_rr_prs_gcd([-2], [2], ZZ) == ([2], [-1], [1])

    assert dup_zz_heu_gcd([2], [-2], ZZ) == ([2], [1], [-1])
    assert dup_rr_prs_gcd([2], [-2], ZZ) == ([2], [1], [-1])

    assert dup_zz_heu_gcd([-2], [-2], ZZ) == ([2], [-1], [-1])
    assert dup_rr_prs_gcd([-2], [-2], ZZ) == ([2], [-1], [-1])

    assert dup_zz_heu_gcd([1,2,1], [1], ZZ) == ([1], [1, 2, 1], [1])
    assert dup_rr_prs_gcd([1,2,1], [1], ZZ) == ([1], [1, 2, 1], [1])

    assert dup_zz_heu_gcd([1,2,1], [2], ZZ) == ([1], [1, 2, 1], [2])
    assert dup_rr_prs_gcd([1,2,1], [2], ZZ) == ([1], [1, 2, 1], [2])

    assert dup_zz_heu_gcd([2,4,2], [2], ZZ) == ([2], [1, 2, 1], [1])
    assert dup_rr_prs_gcd([2,4,2], [2], ZZ) == ([2], [1, 2, 1], [1])

    assert dup_zz_heu_gcd([2], [2,4,2], ZZ) == ([2], [1], [1, 2, 1])
    assert dup_rr_prs_gcd([2], [2,4,2], ZZ) == ([2], [1], [1, 2, 1])

    assert dup_zz_heu_gcd([2,4,2], [1,1], ZZ) == ([1, 1], [2, 2], [1])
    assert dup_rr_prs_gcd([2,4,2], [1,1], ZZ) == ([1, 1], [2, 2], [1])

    assert dup_zz_heu_gcd([1,1], [2,4,2], ZZ) == ([1, 1], [1], [2, 2])
    assert dup_rr_prs_gcd([1,1], [2,4,2], ZZ) == ([1, 1], [1], [2, 2])

    f, g = [1, -31], [1, 0]

    assert dup_zz_heu_gcd(f, g, ZZ) == ([1], f, g)
    assert dup_rr_prs_gcd(f, g, ZZ) == ([1], f, g)

    f = [1,8,21,22,8]
    g = [1,6,11,6]

    h = [1,3,2]

    cff = [1,5,4]
    cfg = [1,3]

    assert dup_zz_heu_gcd(f, g, ZZ) == (h, cff, cfg)
    assert dup_rr_prs_gcd(f, g, ZZ) == (h, cff, cfg)

    f = [1,0,0,0,-4]
    g = [1,0,4,0, 4]

    h = [1,0,2]

    cff = [1,0,-2]
    cfg = [1,0, 2]

    assert dup_zz_heu_gcd(f, g, ZZ) == (h, cff, cfg)
    assert dup_rr_prs_gcd(f, g, ZZ) == (h, cff, cfg)

    f = [1,0,1,0,-3,-3,8,2,-5]
    g = [3,0,5,-0,-4,-9,21]

    h = [1]

    cff = f
    cfg = g

    assert dup_zz_heu_gcd(f, g, ZZ) == (h, cff, cfg)
    assert dup_rr_prs_gcd(f, g, ZZ) == (h, cff, cfg)

    f = dup_normal([1,0,1,0,-3,-3,8,2,-5], QQ)
    g = dup_normal([3,0,5,-0,-4,-9,21], QQ)

    h = dup_normal([1], QQ)

    assert dup_qq_heu_gcd(f, g, QQ) == (h, cff, cfg)
    assert dup_ff_prs_gcd(f, g, QQ) == (h, cff, cfg)

    f = [-352518131239247345597970242177235495263669787845475025293906825864749649589178600387510272,
         0, 0, 0, 0, 0, 0,
         46818041807522713962450042363465092040687472354933295397472942006618953623327997952,
         0, 0, 0, 0, 0, 0,
         378182690892293941192071663536490788434899030680411695933646320291525827756032,
         0, 0, 0, 0, 0, 0,
         112806468807371824947796775491032386836656074179286744191026149539708928,
         0, 0, 0, 0, 0, 0,
         -12278371209708240950316872681744825481125965781519138077173235712,
         0, 0, 0, 0, 0, 0,
         289127344604779611146960547954288113529690984687482920704,
         0, 0, 0, 0, 0, 0,
         19007977035740498977629742919480623972236450681,
         0, 0, 0, 0, 0, 0,
         311973482284542371301330321821976049]

    g = [365431878023781158602430064717380211405897160759702125019136,
         0, 0, 0, 0, 0, 0,
         197599133478719444145775798221171663643171734081650688,
         0, 0, 0, 0, 0, 0,
         -9504116979659010018253915765478924103928886144,
         0, 0, 0, 0, 0, 0,
         -311973482284542371301330321821976049]

    f = dup_normal(f, ZZ)
    g = dup_normal(g, ZZ)

    assert dup_zz_heu_gcd(f, dup_diff(f, 1, ZZ), ZZ)[0] == g
    assert dup_rr_prs_gcd(f, dup_diff(f, 1, ZZ), ZZ)[0] == g

    f = [QQ(1,2),QQ(1),QQ(1,2)]
    g = [QQ(1,2),QQ(1,2)]

    h = [QQ(1), QQ(1)]

    assert dup_qq_heu_gcd(f, g, QQ) == (h, g, [QQ(1,2)])
    assert dup_ff_prs_gcd(f, g, QQ) == (h, g, [QQ(1,2)])

def test_dmp_gcd():
    assert dmp_zz_heu_gcd([[]], [[]], 1, ZZ) == ([[]], [[]], [[]])
    assert dmp_rr_prs_gcd([[]], [[]], 1, ZZ) == ([[]], [[]], [[]])

    assert dmp_zz_heu_gcd([[2]], [[]], 1, ZZ) == ([[2]], [[1]], [[]])
    assert dmp_rr_prs_gcd([[2]], [[]], 1, ZZ) == ([[2]], [[1]], [[]])

    assert dmp_zz_heu_gcd([[-2]], [[]], 1, ZZ) == ([[2]], [[-1]], [[]])
    assert dmp_rr_prs_gcd([[-2]], [[]], 1, ZZ) == ([[2]], [[-1]], [[]])

    assert dmp_zz_heu_gcd([[]], [[-2]], 1, ZZ) == ([[2]], [[]], [[-1]])
    assert dmp_rr_prs_gcd([[]], [[-2]], 1, ZZ) == ([[2]], [[]], [[-1]])

    assert dmp_zz_heu_gcd([[]], [[2],[4]], 1, ZZ) == ([[2],[4]], [[]], [[1]])
    assert dmp_rr_prs_gcd([[]], [[2],[4]], 1, ZZ) == ([[2],[4]], [[]], [[1]])

    assert dmp_zz_heu_gcd([[2],[4]], [[]], 1, ZZ) == ([[2],[4]], [[1]], [[]])
    assert dmp_rr_prs_gcd([[2],[4]], [[]], 1, ZZ) == ([[2],[4]], [[1]], [[]])

    assert dmp_zz_heu_gcd([[2]], [[2]], 1, ZZ) == ([[2]], [[1]], [[1]])
    assert dmp_rr_prs_gcd([[2]], [[2]], 1, ZZ) == ([[2]], [[1]], [[1]])

    assert dmp_zz_heu_gcd([[-2]], [[2]], 1, ZZ) == ([[2]], [[-1]], [[1]])
    assert dmp_rr_prs_gcd([[-2]], [[2]], 1, ZZ) == ([[2]], [[-1]], [[1]])

    assert dmp_zz_heu_gcd([[2]], [[-2]], 1, ZZ) == ([[2]], [[1]], [[-1]])
    assert dmp_rr_prs_gcd([[2]], [[-2]], 1, ZZ) == ([[2]], [[1]], [[-1]])

    assert dmp_zz_heu_gcd([[-2]], [[-2]], 1, ZZ) == ([[2]], [[-1]], [[-1]])
    assert dmp_rr_prs_gcd([[-2]], [[-2]], 1, ZZ) == ([[2]], [[-1]], [[-1]])

    assert dmp_zz_heu_gcd([[1],[2],[1]], [[1]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[1]])
    assert dmp_rr_prs_gcd([[1],[2],[1]], [[1]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[1]])

    assert dmp_zz_heu_gcd([[1],[2],[1]], [[2]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[2]])
    assert dmp_rr_prs_gcd([[1],[2],[1]], [[2]], 1, ZZ) == ([[1]], [[1], [2], [1]], [[2]])

    assert dmp_zz_heu_gcd([[2],[4],[2]], [[2]], 1, ZZ) == ([[2]], [[1], [2], [1]], [[1]])
    assert dmp_rr_prs_gcd([[2],[4],[2]], [[2]], 1, ZZ) == ([[2]], [[1], [2], [1]], [[1]])

    assert dmp_zz_heu_gcd([[2]], [[2],[4],[2]], 1, ZZ) == ([[2]], [[1]], [[1], [2], [1]])
    assert dmp_rr_prs_gcd([[2]], [[2],[4],[2]], 1, ZZ) == ([[2]], [[1]], [[1], [2], [1]])

    assert dmp_zz_heu_gcd([[2],[4],[2]], [[1],[1]], 1, ZZ) == ([[1], [1]], [[2], [2]], [[1]])
    assert dmp_rr_prs_gcd([[2],[4],[2]], [[1],[1]], 1, ZZ) == ([[1], [1]], [[2], [2]], [[1]])

    assert dmp_zz_heu_gcd([[1],[1]], [[2],[4],[2]], 1, ZZ) == ([[1], [1]], [[1]], [[2], [2]])
    assert dmp_rr_prs_gcd([[1],[1]], [[2],[4],[2]], 1, ZZ) == ([[1], [1]], [[1]], [[2], [2]])

    assert dmp_zz_heu_gcd([[[[1,2,1]]]], [[[[2,2]]]], 3, ZZ) == ([[[[1,1]]]], [[[[1,1]]]], [[[[2]]]])
    assert dmp_rr_prs_gcd([[[[1,2,1]]]], [[[[2,2]]]], 3, ZZ) == ([[[[1,1]]]], [[[[1,1]]]], [[[[2]]]])

    f, g = [[[[1,2,1],[1,1],[]]]], [[[[1,2,1]]]]
    h, cff, cfg = [[[[1,1]]]], [[[[1,1],[1],[]]]], [[[[1,1]]]]

    assert dmp_zz_heu_gcd(f, g, 3, ZZ) == (h, cff, cfg)
    assert dmp_rr_prs_gcd(f, g, 3, ZZ) == (h, cff, cfg)

    assert dmp_zz_heu_gcd(g, f, 3, ZZ) == (h, cfg, cff)
    assert dmp_rr_prs_gcd(g, f, 3, ZZ) == (h, cfg, cff)

    f, g, h = dmp_fateman_poly_F_1(2, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    H, cff, cfg = dmp_rr_prs_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    f, g, h = dmp_fateman_poly_F_1(4, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 4, ZZ)

    assert H == h and dmp_mul(H, cff, 4, ZZ) == f \
                  and dmp_mul(H, cfg, 4, ZZ) == g

    f, g, h = dmp_fateman_poly_F_1(6, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 6, ZZ)

    assert H == h and dmp_mul(H, cff, 6, ZZ) == f \
                  and dmp_mul(H, cfg, 6, ZZ) == g

    f, g, h = dmp_fateman_poly_F_1(8, ZZ)

    H, cff, cfg = dmp_zz_heu_gcd(f, g, 8, ZZ)

    assert H == h and dmp_mul(H, cff, 8, ZZ) == f \
                  and dmp_mul(H, cfg, 8, ZZ) == g

    f, g, h = dmp_fateman_poly_F_2(2, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    H, cff, cfg = dmp_rr_prs_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    f, g, h = dmp_fateman_poly_F_3(2, ZZ)
    H, cff, cfg = dmp_zz_heu_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    H, cff, cfg = dmp_rr_prs_gcd(f, g, 2, ZZ)

    assert H == h and dmp_mul(H, cff, 2, ZZ) == f \
                  and dmp_mul(H, cfg, 2, ZZ) == g

    f, g, h = dmp_fateman_poly_F_3(4, ZZ)
    H, cff, cfg = dmp_inner_gcd(f, g, 4, ZZ)

    assert H == h and dmp_mul(H, cff, 4, ZZ) == f \
                  and dmp_mul(H, cfg, 4, ZZ) == g

    f = [[QQ(1,2)],[QQ(1)],[QQ(1,2)]]
    g = [[QQ(1,2)],[QQ(1,2)]]

    h = [[QQ(1)],[QQ(1)]]

    assert dmp_qq_heu_gcd(f, g, 1, QQ) == (h, g, [[QQ(1,2)]])
    assert dmp_ff_prs_gcd(f, g, 1, QQ) == (h, g, [[QQ(1,2)]])

def test_dup_lcm():
    assert dup_lcm([2], [6], ZZ) == [6]

    assert dup_lcm([2,0,0,0], [6,0], ZZ) == [6,0,0,0]
    assert dup_lcm([2,0,0,0], [3,0], ZZ) == [6,0,0,0]

    assert dup_lcm([1,1,0], [1,0], ZZ) == [1,1,0]
    assert dup_lcm([1,1,0], [2,0], ZZ) == [2,2,0]
    assert dup_lcm([1,2,0], [1,0], ZZ) == [1,2,0]
    assert dup_lcm([2,1,0], [1,0], ZZ) == [2,1,0]
    assert dup_lcm([2,1,0], [2,0], ZZ) == [4,2,0]

def test_dmp_lcm():
    assert dmp_lcm([[2]], [[6]], 1, ZZ) == [[6]]
    assert dmp_lcm([[1],[]], [[1,0]], 1, ZZ) == [[1,0],[]]

    assert dmp_lcm([[2],[],[],[]], [[6,0,0],[]], 1, ZZ) == [[6,0,0],[],[],[]]
    assert dmp_lcm([[2],[],[],[]], [[3,0,0],[]], 1, ZZ) == [[6,0,0],[],[],[]]

    assert dmp_lcm([[1,0],[],[]], [[1,0,0],[]], 1, ZZ) == [[1,0,0],[],[]]

    f = [[2,-3,-2,3,0,0],[]]
    g = [[1,0,-2,0,1,0]]
    h = [[2,-3,-4,6,2,-3,0,0],[]]

    assert dmp_lcm(f, g, 1, ZZ) == h

    f = [[1],[-3,0],[-9,0,0],[-5,0,0,0]]
    g = [[1],[6,0],[12,0,0],[10,0,0,0],[3,0,0,0,0]]
    h = [[1],[1,0],[-18,0,0],[-50,0,0,0],[-47,0,0,0,0],[-15,0,0,0,0,0]]

    assert dmp_lcm(f, g, 1, ZZ) == h

def test_dup_trunc():
    assert dup_trunc([1,2,3,4,5,6], ZZ(3), ZZ) == [1, -1, 0, 1, -1, 0]
    assert dup_trunc([6,5,4,3,2,1], ZZ(3), ZZ) ==    [-1, 1, 0, -1, 1]

def test_dmp_trunc():
    assert dmp_trunc([[]], [1,2], 2, ZZ) == [[]]
    assert dmp_trunc([[1,2], [1,4,1], [1]], [1,2], 1, ZZ) == [[-3], [1]]

def test_dmp_ground_trunc():
    assert dmp_ground_trunc(f_0, ZZ(3), 2, ZZ) == \
        dmp_normal([[[1, -1, 0], [-1]], [[]], [[1, -1, 0], [1, -1, 1], [1]]], 2, ZZ)

def test_dup_monic():
    assert dup_monic([3,6,9], ZZ) == [1,2,3]

    raises(ExactQuotientFailed, 'dup_monic([3,4,5], ZZ)')

    assert dup_monic([], QQ) == []
    assert dup_monic([QQ(1)], QQ) == [QQ(1)]
    assert dup_monic([QQ(7),QQ(1),QQ(21)], QQ) == [QQ(1),QQ(1,7),QQ(3)]

def test_dmp_ground_monic():
    assert dmp_ground_monic([[3],[6],[9]], 1, ZZ) == [[1],[2],[3]]

    raises(ExactQuotientFailed, 'dmp_ground_monic([[3],[4],[5]], 1, ZZ)')

    assert dmp_ground_monic([[]], 1, QQ) == [[]]
    assert dmp_ground_monic([[QQ(1)]], 1, QQ) == [[QQ(1)]]
    assert dmp_ground_monic([[QQ(7)],[QQ(1)],[QQ(21)]], 1, QQ) == [[QQ(1)],[QQ(1,7)],[QQ(3)]]

def test_dup_content():
    assert dup_content([], ZZ) == ZZ(0)
    assert dup_content([1], ZZ) == ZZ(1)
    assert dup_content([-1], ZZ) == ZZ(1)
    assert dup_content([1,1], ZZ) == ZZ(1)
    assert dup_content([2,2], ZZ) == ZZ(2)
    assert dup_content([1,2,1], ZZ) == ZZ(1)
    assert dup_content([2,4,2], ZZ) == ZZ(2)

    assert dup_content([QQ(2,3),QQ(4,5)], QQ) == QQ(1)

def test_dmp_content():
    assert dmp_content([[-2]], 1, ZZ) == [2]

    f, g, F = [ZZ(3),ZZ(2),ZZ(1)], [ZZ(1)], []

    for i in xrange(0, 5):
        g = dup_mul(g, f, ZZ)
        F.insert(0, g)

    assert dmp_content(F, 1, ZZ) == f

    assert dmp_one_p(dmp_content(f_4, 2, ZZ), 1, ZZ)
    assert dmp_one_p(dmp_content(f_5, 2, ZZ), 1, ZZ)
    assert dmp_one_p(dmp_content(f_6, 3, ZZ), 2, ZZ)

def test_dmp_ground_content():
    assert dmp_ground_content([[]], 1, ZZ) == ZZ(0)
    assert dmp_ground_content([[]], 1, QQ) == QQ(0)
    assert dmp_ground_content([[1]], 1, ZZ) == ZZ(1)
    assert dmp_ground_content([[-1]], 1, ZZ) == ZZ(1)
    assert dmp_ground_content([[1],[1]], 1, ZZ) == ZZ(1)
    assert dmp_ground_content([[2],[2]], 1, ZZ) == ZZ(2)
    assert dmp_ground_content([[1],[2],[1]], 1, ZZ) == ZZ(1)
    assert dmp_ground_content([[2],[4],[2]], 1, ZZ) == ZZ(2)

    assert dmp_ground_content([[QQ(2,3)],[QQ(4,5)]], 1, QQ) == QQ(1)

    assert dmp_ground_content(f_0, 2, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_0, ZZ(2), 2, ZZ), 2, ZZ) == ZZ(2)

    assert dmp_ground_content(f_1, 2, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_1, ZZ(3), 2, ZZ), 2, ZZ) == ZZ(3)

    assert dmp_ground_content(f_2, 2, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_2, ZZ(4), 2, ZZ), 2, ZZ) == ZZ(4)

    assert dmp_ground_content(f_3, 2, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_3, ZZ(5), 2, ZZ), 2, ZZ) == ZZ(5)

    assert dmp_ground_content(f_4, 2, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_4, ZZ(6), 2, ZZ), 2, ZZ) == ZZ(6)

    assert dmp_ground_content(f_5, 2, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_5, ZZ(7), 2, ZZ), 2, ZZ) == ZZ(7)

    assert dmp_ground_content(f_6, 3, ZZ) == ZZ(1)
    assert dmp_ground_content(dmp_mul_ground(f_6, ZZ(8), 3, ZZ), 3, ZZ) == ZZ(8)

def test_dup_primitive():
    assert dup_primitive([], ZZ) == (0, [])
    assert dup_primitive([1], ZZ) == (1, [1])
    assert dup_primitive([1,1], ZZ) == (1, [1,1])
    assert dup_primitive([2,2], ZZ) == (2, [1,1])
    assert dup_primitive([1,2,1], ZZ) == (1, [1,2,1])
    assert dup_primitive([2,4,2], ZZ) == (2, [1,2,1])

    assert dup_primitive([QQ(2,3),QQ(4,5)], QQ) == (QQ(1), [QQ(2,3),QQ(4,5)])

def test_dmp_primitive():
    assert dmp_primitive([[]], 1, ZZ) == ([], [[]])
    assert dmp_primitive([[1]], 1, ZZ) == ([1], [[1]])

    f, g, F = [ZZ(3),ZZ(2),ZZ(1)], [ZZ(1)], []

    for i in xrange(0, 5):
        g = dup_mul(g, f, ZZ)
        F.insert(0, g)

    assert dmp_primitive(F, 1, ZZ) == (f,
        [ dup_exquo(c, f, ZZ) for c in F ])

    cont, f = dmp_primitive(f_4, 2, ZZ)
    assert dmp_one_p(cont, 1, ZZ) and f == f_4
    cont, f = dmp_primitive(f_5, 2, ZZ)
    assert dmp_one_p(cont, 1, ZZ) and f == f_5
    cont, f = dmp_primitive(f_6, 3, ZZ)
    assert dmp_one_p(cont, 2, ZZ) and f == f_6

def test_dmp_ground_primitive():
    assert dmp_ground_primitive([[]], 1, ZZ) == (ZZ(0), [[]])

    assert dmp_ground_primitive(f_0, 2, ZZ) == (ZZ(1), f_0)
    assert dmp_ground_primitive(dmp_mul_ground(f_0, ZZ(2), 2, ZZ), 2, ZZ) == (ZZ(2), f_0)

    assert dmp_ground_primitive(f_1, 2, ZZ) == (ZZ(1), f_1)
    assert dmp_ground_primitive(dmp_mul_ground(f_1, ZZ(3), 2, ZZ), 2, ZZ) == (ZZ(3), f_1)

    assert dmp_ground_primitive(f_2, 2, ZZ) == (ZZ(1), f_2)
    assert dmp_ground_primitive(dmp_mul_ground(f_2, ZZ(4), 2, ZZ), 2, ZZ) == (ZZ(4), f_2)

    assert dmp_ground_primitive(f_3, 2, ZZ) == (ZZ(1), f_3)
    assert dmp_ground_primitive(dmp_mul_ground(f_3, ZZ(5), 2, ZZ), 2, ZZ) == (ZZ(5), f_3)

    assert dmp_ground_primitive(f_4, 2, ZZ) == (ZZ(1), f_4)
    assert dmp_ground_primitive(dmp_mul_ground(f_4, ZZ(6), 2, ZZ), 2, ZZ) == (ZZ(6), f_4)

    assert dmp_ground_primitive(f_5, 2, ZZ) == (ZZ(1), f_5)
    assert dmp_ground_primitive(dmp_mul_ground(f_5, ZZ(7), 2, ZZ), 2, ZZ) == (ZZ(7), f_5)

    assert dmp_ground_primitive(f_6, 3, ZZ) == (ZZ(1), f_6)
    assert dmp_ground_primitive(dmp_mul_ground(f_6, ZZ(8), 3, ZZ), 3, ZZ) == (ZZ(8), f_6)

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

def test_dup_extract():
    f = dup_normal([2930944, 0, 2198208, 0, 549552, 0, 45796], ZZ)
    g = dup_normal([17585664, 0, 8792832, 0, 1099104, 0], ZZ)

    F = dup_normal([64, 0, 48, 0, 12, 0, 1], ZZ)
    G = dup_normal([384, 0, 192, 0, 24, 0], ZZ)

    assert dup_extract(f, g, ZZ) == (45796, F, G)

def test_dmp_ground_extract():
    f = dmp_normal([[2930944], [], [2198208], [], [549552], [], [45796]], 1, ZZ)
    g = dmp_normal([[17585664], [], [8792832], [], [1099104], []], 1, ZZ)

    F = dmp_normal([[64], [], [48], [], [12], [], [1]], 1, ZZ)
    G = dmp_normal([[384], [], [192], [], [24], []], 1, ZZ)

    assert dmp_ground_extract(f, g, 1, ZZ) == (45796, F, G)

def test_dup_real_imag():
    assert dup_real_imag([], ZZ) == ([[]], [[]])
    assert dup_real_imag([1], ZZ) == ([[1]], [[]])

    assert dup_real_imag([1,1], ZZ) == ([[1], [1]], [[1,0]])
    assert dup_real_imag([1,2], ZZ) == ([[1], [2]], [[1,0]])

    assert dup_real_imag([1,2,3], ZZ) == ([[1], [2], [-1,0,3]], [[2,0], [2,0]])

    raises(DomainError, "dup_real_imag([EX(1), EX(2)], EX)")

def test_dup_mirror():
    assert dup_mirror([], ZZ) == []
    assert dup_mirror([1], ZZ) == [1]

    assert dup_mirror([1,2,3,4,5], ZZ) == [1,-2,3,-4,5]
    assert dup_mirror([1,2,3,4,5,6], ZZ) == [-1,2,-3,4,-5,6]

def test_dup_scale():
    assert dup_scale([], -1, ZZ) == []
    assert dup_scale([1], -1, ZZ) == [1]

    assert dup_scale([1,2,3,4,5], -1, ZZ) == [1,-2,3,-4,5]
    assert dup_scale([1,2,3,4,5], -7, ZZ) == [2401,-686,147,-28,5]

def test_dup_taylor():
    assert dup_taylor([], 1, ZZ) == []
    assert dup_taylor([1], 1, ZZ) == [1]

    assert dup_taylor([1,2,3,4,5], 1, ZZ) == [1,6,15,20,15]
    assert dup_taylor([1,2,3,4,5], 7, ZZ) == [1,30,339,1712,3267]

def test_dup_transform():
    assert dup_transform([], [], [1,1], ZZ) == []
    assert dup_transform([], [1], [1,1], ZZ) == []
    assert dup_transform([], [1,2], [1,1], ZZ) == []

    assert dup_transform([6,-5,4,-3,17], [1,-3,4], [2,-3], ZZ) == \
        [6,-82,541,-2205,6277,-12723,17191,-13603,4773]

def test_dup_compose():
    assert dup_compose([], [], ZZ) == []
    assert dup_compose([], [1], ZZ) == []
    assert dup_compose([], [1,2], ZZ) == []

    assert dup_compose([1], [], ZZ) == [1]

    assert dup_compose([1,2,0], [], ZZ) == []
    assert dup_compose([1,2,1], [], ZZ) == [1]

    assert dup_compose([1,2,1], [1], ZZ) == [4]
    assert dup_compose([1,2,1], [7], ZZ) == [64]

    assert dup_compose([1,2,1], [1,-1], ZZ) == [1,0,0]
    assert dup_compose([1,2,1], [1, 1], ZZ) == [1,4,4]
    assert dup_compose([1,2,1], [1, 2, 1], ZZ) == [1,4,8,8,4]

def test_dmp_compose():
    assert dmp_compose([1,2,1], [1,2,1], 0, ZZ) == [1,4,8,8,4]

    assert dmp_compose([[[]]], [[[]]], 2, ZZ) == [[[]]]
    assert dmp_compose([[[]]], [[[1]]], 2, ZZ) == [[[]]]
    assert dmp_compose([[[]]], [[[1]],[[2]]], 2, ZZ) == [[[]]]

    assert dmp_compose([[[1]]], [], 2, ZZ) == [[[1]]]

    assert dmp_compose([[1],[2],[ ]], [[]], 1, ZZ) == [[]]
    assert dmp_compose([[1],[2],[1]], [[]], 1, ZZ) == [[1]]

    assert dmp_compose([[1],[2],[1]], [[1]], 1, ZZ) == [[4]]
    assert dmp_compose([[1],[2],[1]], [[7]], 1, ZZ) == [[64]]

    assert dmp_compose([[1],[2],[1]], [[1],[-1]], 1, ZZ) == [[1],[ ],[ ]]
    assert dmp_compose([[1],[2],[1]], [[1],[ 1]], 1, ZZ) == [[1],[4],[4]]

    assert dmp_compose([[1],[2],[1]], [[1], [2], [1]], 1, ZZ) == [[1],[4],[8],[8],[4]]

def test_dup_decompose():
    assert dup_decompose([1], ZZ) == ([1],)

    assert dup_decompose([1,0], ZZ) == ([1,0],)
    assert dup_decompose([1,0,0,0], ZZ) == ([1,0,0,0],)

    assert dup_decompose([1,0,0,0,0], ZZ) == ([1,0,0], [1,0,0])
    assert dup_decompose([1,0,0,0,0,0,0], ZZ) == ([1,0,0,0], [1,0,0])

    assert dup_decompose([7,0,0,0,1], ZZ) == ([7,0,1], [1,0,0])
    assert dup_decompose([4,0,3,0,2], ZZ) == ([4,3,2], [1,0,0])

    f = [1,0,20,0,150,0,500,0,625,-2,0,-10,9]

    assert dup_decompose(f, ZZ) == ([1,0,0,-2,9], [1,0,5,0])

    f = [2,0,40,0,300,0,1000,0,1250,-4,0,-20,18]

    assert dup_decompose(f, ZZ) == ([2,0,0,-4,18], [1,0,5,0])

    f = [1,0,20,-8,150,-120,524,-600,865,-1034,600,-170,29]

    assert dup_decompose(f, ZZ) == ([1,-8,24,-34,29], [1,0,5,0])

    f = [DMP([6,0,-42], ZZ), DMP([48,0,96], ZZ), DMP([144,648,288], ZZ),
         DMP([624,864,384], ZZ), DMP([108,312,432,192], ZZ)]

    assert dup_decompose(f, ZZ['a']) == (f,)

def test_dup_sturm():
    assert dup_sturm([QQ(5)], QQ) == [[QQ(1)]]
    assert dup_sturm([QQ(1),QQ(0)], QQ) == [[QQ(1),QQ(0)], [QQ(1)]]

    f = dup_normal([1,-2,3,-5], QQ)

    assert dup_sturm(f, QQ) == \
        [f, [QQ(3),QQ(-4),QQ(3)], [QQ(-10,9),QQ(13,3)], [QQ(-3303,100)]]

def test_dmp_lift():
    q = [QQ(1,1),QQ(0,1),QQ(1,1)]

    f = [ANP([QQ(1,1)], q, QQ), ANP([], q, QQ), ANP([], q, QQ),
         ANP([QQ(1,1),QQ(0,1)], q, QQ), ANP([QQ(17,1),QQ(0,1)], q, QQ)]

    assert dmp_lift(f, 0, QQ.algebraic_field(I)) == \
        [QQ(1),QQ(0),QQ(0),QQ(0),QQ(0),QQ(0),QQ(2),QQ(0),QQ(578),
         QQ(0),QQ(0),QQ(0),QQ(1),QQ(0),QQ(-578),QQ(0),QQ(83521)]

    raises(DomainError, "dmp_lift([EX(1), EX(2)], 0, EX)")

def test_dup_sign_variations():
    assert dup_sign_variations([], ZZ) == 0
    assert dup_sign_variations([1,0], ZZ) == 0
    assert dup_sign_variations([1,0,2], ZZ) == 0
    assert dup_sign_variations([1,0,3,0], ZZ) == 0
    assert dup_sign_variations([1,0,4,0,5], ZZ) == 0

    assert dup_sign_variations([-1,0,2], ZZ) == 1
    assert dup_sign_variations([-1,0,3,0], ZZ) == 1
    assert dup_sign_variations([-1,0,4,0,5], ZZ) == 1

    assert dup_sign_variations([-1,-4,-5], ZZ) == 0
    assert dup_sign_variations([ 1,-4,-5], ZZ) == 1
    assert dup_sign_variations([ 1, 4,-5], ZZ) == 1
    assert dup_sign_variations([ 1,-4, 5], ZZ) == 2
    assert dup_sign_variations([-1, 4,-5], ZZ) == 2
    assert dup_sign_variations([-1, 4, 5], ZZ) == 1
    assert dup_sign_variations([-1,-4, 5], ZZ) == 1
    assert dup_sign_variations([ 1, 4, 5], ZZ) == 0

    assert dup_sign_variations([-1,0,-4,0,-5], ZZ) == 0
    assert dup_sign_variations([ 1,0,-4,0,-5], ZZ) == 1
    assert dup_sign_variations([ 1,0, 4,0,-5], ZZ) == 1
    assert dup_sign_variations([ 1,0,-4,0, 5], ZZ) == 2
    assert dup_sign_variations([-1,0, 4,0,-5], ZZ) == 2
    assert dup_sign_variations([-1,0, 4,0, 5], ZZ) == 1
    assert dup_sign_variations([-1,0,-4,0, 5], ZZ) == 1
    assert dup_sign_variations([ 1,0, 4,0, 5], ZZ) == 0

