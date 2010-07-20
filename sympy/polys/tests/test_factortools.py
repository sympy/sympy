"""Tools for polynomial factorization routines in characteristic zero. """

from sympy.polys.densearith import (
    dup_mul_ground, dmp_mul_ground,
    dup_pow, dmp_pow,
    dmp_expand,
)

from sympy.polys.densebasic import (
    dup_degree, dmp_degree,
    dup_normal, dmp_normal,
    dup_from_raw_dict, dup_to_raw_dict,
    dmp_from_dict, dmp_to_dict,
    dmp_nest, dmp_raise,
)

from sympy.polys.densetools import (
    dup_primitive, dup_sqf_p, dmp_eval_tail,
)

from sympy.polys.factortools import (
    dup_trial_division, dmp_trial_division,
    dup_zz_mignotte_bound, dmp_zz_mignotte_bound,
    dup_zz_hensel_step, dup_zz_hensel_lift,
    dup_zz_irreducible_p,
    dup_zz_zassenhaus, dmp_zz_wang,
    dmp_zz_wang_non_divisors,
    dmp_zz_wang_test_points,
    dmp_zz_wang_lead_coeffs,
    dmp_zz_wang_hensel_lifting,
    dup_zz_diophantine, dmp_zz_diophantine,
    dup_zz_cyclotomic_poly, dup_zz_cyclotomic_factor,
    dup_zz_factor, dup_zz_factor_sqf, dmp_zz_factor,
    dup_ext_factor, dmp_ext_factor,
    dup_factor_list, dmp_factor_list,
    dup_factor_list_include, dmp_factor_list_include,
)

from sympy.polys.specialpolys import (
    f_1, f_2, f_3, f_4, f_5, f_6, w_1, w_2,
)

from sympy.polys.polyconfig import setup
from sympy.polys.polyerrors import DomainError
from sympy.polys.polyclasses import DMP, DMF, ANP
from sympy.polys.domains import FF, ZZ, QQ, RR, EX

from sympy import raises, nextprime, sin, sqrt, I

def test_dup_trial_division():
    assert dup_trial_division([1,8,25,38,28,8],
        ([1,1], [1,2]), ZZ) == [([1,1], 2), ([1,2], 3)]

def test_dmp_trial_division():
    assert dmp_trial_division([[1],[8],[25],[38],[28],[8]],
        ([[1],[1]], [[1],[2]]), 1, ZZ) == [([[1],[1]], 2), ([[1],[2]], 3)]

def test_dup_zz_mignotte_bound():
    assert dup_zz_mignotte_bound([2,3,4], ZZ) == 32

def test_dmp_zz_mignotte_bound():
    assert dmp_zz_mignotte_bound([[2],[3],[4]], 1, ZZ) == 32

def test_dup_zz_hensel_step():
    f = dup_from_raw_dict({4:1, 0:-1}, ZZ)

    g = dup_from_raw_dict({3:1, 2:2, 1:-1, 0:-2}, ZZ)
    h = dup_from_raw_dict({1:1, 0:-2}, ZZ)
    s = dup_from_raw_dict({0:-2}, ZZ)
    t = dup_from_raw_dict({2:2, 1:-2, 0:-1}, ZZ)

    G, H, S, T = dup_zz_hensel_step(5, f, g, h, s, t, ZZ)

    assert G == dup_from_raw_dict({3:1, 2:7, 1:-1, 0:-7}, ZZ)
    assert H == dup_from_raw_dict({1:1, 0:-7}, ZZ)
    assert S == dup_from_raw_dict({0:8}, ZZ)
    assert T == dup_from_raw_dict({2:-8, 1:-12, 0:-1}, ZZ)

def test_dup_zz_hensel_lift():
    f = dup_from_raw_dict({4:1, 0:-1}, ZZ)

    f1 = dup_from_raw_dict({1:1, 0:-1}, ZZ)
    f2 = dup_from_raw_dict({1:1, 0:-2}, ZZ)
    f3 = dup_from_raw_dict({1:1, 0: 2}, ZZ)
    f4 = dup_from_raw_dict({1:1, 0: 1}, ZZ)

    ff_list = dup_zz_hensel_lift(5, f, [f1, f2, f3, f4], 4, ZZ)

    assert dup_to_raw_dict(ff_list[0]) == {0: -1,   1: 1}
    assert dup_to_raw_dict(ff_list[1]) == {0: -182, 1: 1}
    assert dup_to_raw_dict(ff_list[2]) == {0:  182, 1: 1}
    assert dup_to_raw_dict(ff_list[3]) == {0:  1,   1: 1}

def test_dup_zz_irreducible_p():
    assert dup_zz_irreducible_p([3, 2, 6, 8, 7], ZZ) is None
    assert dup_zz_irreducible_p([3, 2, 6, 8, 4], ZZ) is None

    assert dup_zz_irreducible_p([3, 2, 6, 8, 10], ZZ) == True
    assert dup_zz_irreducible_p([3, 2, 6, 8, 14], ZZ) == True

def test_dup_zz_cyclotomic_poly():
    assert dup_zz_cyclotomic_poly(1, ZZ) == [1,-1]
    assert dup_zz_cyclotomic_poly(2, ZZ) == [1,1]
    assert dup_zz_cyclotomic_poly(3, ZZ) == [1,1,1]
    assert dup_zz_cyclotomic_poly(4, ZZ) == [1,0,1]
    assert dup_zz_cyclotomic_poly(5, ZZ) == [1,1,1,1,1]
    assert dup_zz_cyclotomic_poly(6, ZZ) == [1,-1,1]
    assert dup_zz_cyclotomic_poly(7, ZZ) == [1,1,1,1,1,1,1]
    assert dup_zz_cyclotomic_poly(8, ZZ) == [1,0,0,0,1]
    assert dup_zz_cyclotomic_poly(9, ZZ) == [1,0,0,1,0,0,1]

def test_dup_zz_cyclotomic_factor():
    assert dup_zz_cyclotomic_factor([], ZZ) is None
    assert dup_zz_cyclotomic_factor([1], ZZ) is None

    f = dup_from_raw_dict({10:2, 0:-1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) is None
    f = dup_from_raw_dict({10:1, 0:-3}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) is None
    f = dup_from_raw_dict({10:1, 5:1, 0:-1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) is None

    f = dup_from_raw_dict({1:1,0:1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) == \
         [[1, 1]]

    f = dup_from_raw_dict({1:1,0:-1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) == \
        [[1, -1]]

    f = dup_from_raw_dict({2:1,0:1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) == \
        [[1, 0, 1]]

    f = dup_from_raw_dict({2:1,0:-1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) == \
        [[1,-1],
         [1, 1]]

    f = dup_from_raw_dict({27:1,0:1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) == \
        [[1, 1],
         [1, -1, 1],
         [1, 0, 0, -1, 0, 0, 1],
         [1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1]]

    f = dup_from_raw_dict({27:1,0:-1}, ZZ)
    assert dup_zz_cyclotomic_factor(f, ZZ) == \
        [[1, -1],
         [1, 1, 1],
         [1, 0, 0, 1, 0, 0, 1],
         [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]]

def test_dup_zz_factor():
    assert dup_zz_factor([], ZZ) == (0, [])
    assert dup_zz_factor([7], ZZ) == (7, [])
    assert dup_zz_factor([-7], ZZ) == (-7, [])

    assert dup_zz_factor_sqf([], ZZ) == (0, [])
    assert dup_zz_factor_sqf([7], ZZ) == (7, [])
    assert dup_zz_factor_sqf([-7], ZZ) == (-7, [])

    assert dup_zz_factor([2,4], ZZ) == \
        (2, [([1, 2], 1)])
    assert dup_zz_factor_sqf([2,4], ZZ) == \
        (2, [([1, 2], 1)])

    f = [1,0,0,1,1]

    for i in xrange(0, 20):
        assert dup_zz_factor(f, ZZ) == (1, [(f, 1)])

    assert dup_zz_factor([1,2,2], ZZ) == \
        (1, [([1,2,2], 1)])

    assert dup_zz_factor([18,12,2], ZZ) == \
        (2, [([3, 1], 2)])

    assert dup_zz_factor([-9,0,1], ZZ) == \
        (-1, [([3,-1], 1),
              ([3, 1], 1)])

    assert dup_zz_factor_sqf([-9,0,1], ZZ) == \
        (-1, [[3,-1],
              [3, 1]])

    assert dup_zz_factor([1,-6,11,-6], ZZ) == \
        (1, [([1,-3], 1),
             ([1,-2], 1),
             ([1,-1], 1)])

    assert dup_zz_factor_sqf([1,-6,11,-6], ZZ) == \
        (1, [[1,-3],
             [1,-2],
             [1,-1]])

    assert dup_zz_factor([3,10,13,10], ZZ) == \
        (1, [([1,2], 1),
             ([3,4,5], 1)])

    assert dup_zz_factor_sqf([3,10,13,10], ZZ) == \
        (1, [[1,2],
             [3,4,5]])

    assert dup_zz_factor([-1,0,0,0,1,0,0], ZZ) == \
        (-1, [([1,-1], 1),
              ([1, 1], 1),
              ([1, 0], 2),
              ([1, 0, 1], 1)])

    f = [1080, 5184, 2099, 744, 2736, -648, 129, 0, -324]

    assert dup_zz_factor(f, ZZ) == \
        (1, [([5, 24, 9, 0, 12], 1),
             ([216, 0, 31, 0, -27], 1)])

    f = [-29802322387695312500000000000000000000,
          0, 0, 0, 0,
          2980232238769531250000000000000000,
          0, 0, 0, 0,
          1743435859680175781250000000000,
          0, 0, 0, 0,
          114142894744873046875000000,
          0, 0, 0, 0,
          -210106372833251953125,
          0, 0, 0, 0,
          95367431640625]

    assert dup_zz_factor(f, ZZ) == \
        (-95367431640625, [([5, -1], 1),
                           ([100, 10, -1], 2),
                           ([625, 125, 25, 5, 1], 1),
                           ([10000, -3000, 400, -20, 1], 2),
                           ([10000,  2000, 400,  30, 1], 2)])

    f = dup_from_raw_dict({10:1, 0:-1}, ZZ)

    setup('USE_CYCLOTOMIC_FACTOR', True)
    F_0 = dup_zz_factor(f, ZZ)

    setup('USE_CYCLOTOMIC_FACTOR', False)
    F_1 = dup_zz_factor(f, ZZ)

    assert F_0 == F_1 == \
        (1, [([1,-1], 1),
             ([1, 1], 1),
             ([1,-1, 1,-1, 1], 1),
             ([1, 1, 1, 1, 1], 1)])

    setup('USE_CYCLOTOMIC_FACTOR')

    f = dup_from_raw_dict({10:1, 0:1}, ZZ)

    setup('USE_CYCLOTOMIC_FACTOR', True)
    F_0 = dup_zz_factor(f, ZZ)

    setup('USE_CYCLOTOMIC_FACTOR', False)
    F_1 = dup_zz_factor(f, ZZ)

    assert F_0 == F_1 == \
        (1, [([1, 0, 1], 1),
             ([1, 0, -1, 0, 1, 0, -1, 0, 1], 1)])

    setup('USE_CYCLOTOMIC_FACTOR')

def test_dmp_zz_wang():
    p = ZZ(nextprime(dmp_zz_mignotte_bound(w_1, 2, ZZ)))

    assert p == ZZ(6291469)

    t_1, k_1, e_1 = dmp_normal([[1],[]], 1, ZZ), 1, ZZ(-14)
    t_2, k_2, e_2 = dmp_normal([[1, 0]], 1, ZZ), 2, ZZ(3)
    t_3, k_3, e_3 = dmp_normal([[1],[ 1, 0]], 1, ZZ), 2, ZZ(-11)
    t_4, k_4, e_4 = dmp_normal([[1],[-1, 0]], 1, ZZ), 1, ZZ(-17)

    T = [t_1, t_2, t_3, t_4]
    K = [k_1, k_2, k_3, k_4]
    E = [e_1, e_2, e_3, e_4]

    T = zip(T, K)

    A = [ZZ(-14), ZZ(3)]

    S = dmp_eval_tail(w_1, A, 2, ZZ)
    cs, s = dup_primitive(S, ZZ)

    assert cs == 1 and s == S == \
        dup_normal([1036728, 915552, 55748, 105621, -17304, -26841, -644], ZZ)

    assert dmp_zz_wang_non_divisors(E, cs, 4, ZZ) == [7, 3, 11, 17]
    assert dup_sqf_p(s, ZZ) and dup_degree(s) == dmp_degree(w_1, 2)

    _, H = dup_zz_factor_sqf(s, ZZ)

    h_1 = dup_normal([44,  42,   1], ZZ)
    h_2 = dup_normal([126, -9,  28], ZZ)
    h_3 = dup_normal([187,  0, -23], ZZ)

    assert H == [h_1, h_2, h_3]

    lc_1 = dmp_normal([[-4], [-4,0]], 1, ZZ)
    lc_2 = dmp_normal([[-1,0,0], []], 1, ZZ)
    lc_3 = dmp_normal([[1], [], [-1,0,0]], 1, ZZ)

    LC = [lc_1, lc_2, lc_3]

    assert dmp_zz_wang_lead_coeffs(w_1, T, cs, E, H, A, 2, ZZ) == (w_1, H, LC)

    H_1 = [ dmp_normal(t, 0, ZZ) for t in [[44L,42L,1L],[126L,-9L,28L],[187L,0L,-23L]] ]
    H_2 = [ dmp_normal(t, 1, ZZ) for t in [[[-4,-12],[-3,0],[1]],[[-9,0],[-9],[-2,0]],[[1,0,-9],[],[1,-9]]] ]
    H_3 = [ dmp_normal(t, 1, ZZ) for t in [[[-4,-12],[-3,0],[1]],[[-9,0],[-9],[-2,0]],[[1,0,-9],[],[1,-9]]] ]

    c_1 = dmp_normal([-70686,-5863,-17826,2009,5031,74], 0, ZZ)
    c_2 = dmp_normal([[9,12,-45,-108,-324],[18,-216,-810,0],[2,9,-252,-288,-945],[-30,-414,0],[2,-54,-3,81],[12,0]], 1, ZZ)
    c_3 = dmp_normal([[-36,-108,0],[-27,-36,-108],[-8,-42,0],[-6,0,9],[2,0]], 1, ZZ)

    T_1 = [ dmp_normal(t, 0, ZZ) for t in [[-3,0],[-2],[1]] ]
    T_2 = [ dmp_normal(t, 1, ZZ) for t in [[[-1,0],[]],[[-3],[]],[[-6]]] ]
    T_3 = [ dmp_normal(t, 1, ZZ) for t in [[[]],[[]],[[-1]]] ]

    assert dmp_zz_diophantine(H_1, c_1,        [], 5, p, 0, ZZ) == T_1
    assert dmp_zz_diophantine(H_2, c_2, [ZZ(-14)], 5, p, 1, ZZ) == T_2
    assert dmp_zz_diophantine(H_3, c_3, [ZZ(-14)], 5, p, 1, ZZ) == T_3

    factors = dmp_zz_wang_hensel_lifting(w_1, H, LC, A, p, 2, ZZ)

    assert dmp_expand(factors, 2, ZZ) == w_1

def test_dmp_zz_factor():
    assert dmp_zz_factor([], 0, ZZ) == (0, [])
    assert dmp_zz_factor([7], 0, ZZ) == (7, [])
    assert dmp_zz_factor([-7], 0, ZZ) == (-7, [])

    assert dmp_zz_factor([[]], 1, ZZ) == (0, [])
    assert dmp_zz_factor([[7]], 1, ZZ) == (7, [])
    assert dmp_zz_factor([[-7]], 1, ZZ) == (-7, [])

    assert dmp_zz_factor([[1], []], 1, ZZ) == \
        (1, [([[1], []], 1)])

    assert dmp_zz_factor([[4], []], 1, ZZ) == \
        (4, [([[1], []], 1)])

    assert dmp_zz_factor([[4], [2]], 1, ZZ) == \
        (2, [([[2], [1]], 1)])

    assert dmp_zz_factor([[1, 0], [1]], 1, ZZ) == \
        (1, [([[1, 0], [1]], 1)])

    assert dmp_zz_factor([[1,0,1]], 1, ZZ) == \
        (1, [([[1, 0, 1]], 1)])

    assert dmp_zz_factor([[1,0,-1]], 1, ZZ) == \
        (1, [([[1,-1]], 1),
             ([[1, 1]], 1)])

    assert dmp_zz_factor([[1, 6, 9], [], [-1]], 1, ZZ) == \
        (1, [([[1, 3], [-1]], 1), ([[1, 3], [1]], 1)])

    assert dmp_zz_factor([1, 0, -9], 0, ZZ) == \
        (1, [([1, -3], 1), ([1, 3], 1)])

    assert dmp_zz_factor([[1, 0, 0], [], [-9]], 1, ZZ) == \
        (1, [([[1, 0], [-3]], 1), ([[1, 0], [3]], 1)])

    assert dmp_zz_factor([[[1, 0, 0], [], []], [[]], [[-9]]], 2, ZZ) == \
        (1, [([[[1, 0], []], [[-3]]], 1), ([[[1, 0], []], [[3]]], 1)])

    assert dmp_zz_factor([[[[1, 0, 0], [], []], [[]], [[]]], [[[]]], [[[-9]]]], 3, ZZ) == \
        (1, [([[[[1, 0], []], [[]]], [[[-3]]]], 1), ([[[[1, 0], []], [[]]], [[[3]]]], 1)])

    assert dmp_zz_factor(f_1, 2, ZZ) == \
        (1, [([[[1]], [[1, 0], [20]]], 1),
             ([[[1], []], [[1, 10]]],  1),
             ([[[1, 0]], [[1], [30]]], 1)])

    assert dmp_zz_factor(f_2, 2, ZZ) == \
        (1, [([[[1], [], [1, 0, 0]], [[]], [[1], [90]]], 1),
             ([[[1], [1, 0]], [[]], [[]], [[1, -11]]],   1)])

    assert dmp_zz_factor(f_3, 2, ZZ) == \
        (1, [([[[1], [], []], [[1, 0, 0, 0, 1]], [[1, 0]]], 1),
             ([[[1]], [[]], [[1, 0], []], [[1], [1, 0, 0, 0], []]], 1)])

    assert dmp_zz_factor(f_4, 2, ZZ) == \
        (-1, [([[[1], [], [], []], [[1, 0, 0]]], 1),
              ([[[1, 0]], [[]], [[1, 0, 0], [], [], [], [5]]], 1),
              ([[[1], []], [[]], [[]], [[-1, 0, -3]]], 1),
              ([[[1], [], [], [], []], [[]], [[]], [[1, 0, 0]]], 1)])

    assert dmp_zz_factor(f_5, 2, ZZ) == \
        (-1, [([[[1]], [[1], [-1, 0]]], 3)])

    assert dmp_zz_factor(f_6, 3, ZZ) == \
        (1, [([[[[47]], [[]]], [[[1, 0, 0], [], [], [-1, 0, 0]]]], 1),
             ([[[[45]]], [[[]]], [[[]]], [[[-9]], [[-1]], [[]], [[3], [], [2, 0], []]]], 1)])

    assert dmp_zz_factor(w_1, 2, ZZ) == \
        (1, [([[[1], [], [-1, 0, 0]], [[]], [[1], [-1, 0, 0]]], 1),
             ([[[1, 0, 0], []], [[3, 0]], [[2], []]], 1),
             ([[[4], [4, 0]], [[1, 0], []], [[-1]]], 1)])

    f = [[-12, 0], [],
         [], [],
         [240, 0, 0, 0], [],
         [-768, 0, 0, 0, 0], [],
         [1080, 0, 0, 0, 0, 0], [],
         [-768, 0, 0, 0, 0, 0, 0], [],
         [240, 0, 0, 0, 0, 0, 0, 0], [],
         [], [],
         [-12, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

    assert dmp_zz_factor(f, 1, ZZ) == \
        (-12, [([[1, 0]], 1),
               ([[1], [], [-1, 0]], 6),
               ([[1], [], [6, 0], [], [1, 0, 0]], 1)])

def test_dup_ext_factor():
    h = [QQ(1),QQ(0),QQ(1)]
    K = QQ.algebraic_field(I)

    assert dup_ext_factor([], K) == (ANP([], h, QQ), [])

    f = [ANP([QQ(1)], h, QQ), ANP([QQ(1)], h, QQ)]

    assert dup_ext_factor(f, K) == (ANP([QQ(1)], h, QQ), [(f, 1)])

    g = [ANP([QQ(2)], h, QQ), ANP([QQ(2)], h, QQ)]

    assert dup_ext_factor(g, K) == (ANP([QQ(2)], h, QQ), [(f, 1)])

    f = [ANP([QQ(7)], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([QQ(1,1)], h, QQ)]
    g = [ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([QQ(1,7)], h, QQ)]

    assert dup_ext_factor(f, K) == (ANP([QQ(7)], h, QQ), [(g, 1)])

    f = [ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([QQ(1)], h, QQ)]

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(1,1)], h, QQ), [
            ([ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([QQ(-1),QQ(0)], h, QQ)], 1),
            ([ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([QQ( 1),QQ(0)], h, QQ)], 1),
         ])

    f = [ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([QQ(1)], h, QQ)]

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(1,1)], h, QQ), [
            ([ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([QQ(-1),QQ(0)], h, QQ)], 1),
            ([ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([QQ( 1),QQ(0)], h, QQ)], 1),
         ])

    h = [QQ(1),QQ(0),QQ(-2)]
    K = QQ.algebraic_field(sqrt(2))

    f = [ANP([QQ(1)], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([], h, QQ), ANP([QQ(1,1)], h, QQ)]

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(1)], h, QQ), [
            ([ANP([QQ(1)], h, QQ), ANP([QQ(-1),QQ(0)], h, QQ), ANP([QQ(1)], h, QQ)], 1),
            ([ANP([QQ(1)], h, QQ), ANP([QQ( 1),QQ(0)], h, QQ), ANP([QQ(1)], h, QQ)], 1),
         ])

    f = [ANP([QQ(1,1)], h, QQ), ANP([2,0], h, QQ), ANP([QQ(2,1)], h, QQ)]

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(1,1)], h, QQ), [
            ([ANP([1], h, QQ), ANP([1,0], h, QQ)], 2),
        ])

    assert dup_ext_factor(dup_pow(f, 3, K), K) == \
        (ANP([QQ(1,1)], h, QQ), [
            ([ANP([1], h, QQ), ANP([1,0], h, QQ)], 6),
        ])

    f = dup_mul_ground(f, ANP([QQ(2,1)], h, QQ), K)

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(2,1)], h, QQ), [
            ([ANP([1], h, QQ), ANP([1,0], h, QQ)], 2),
        ])

    assert dup_ext_factor(dup_pow(f, 3, K), K) == \
        (ANP([QQ(8,1)], h, QQ), [
            ([ANP([1], h, QQ), ANP([1,0], h, QQ)], 6),
        ])

    h = [QQ(1,1), QQ(0,1), QQ(1,1)]
    K = QQ.algebraic_field(I)

    f = [ANP([QQ(4,1)], h, QQ), ANP([], h, QQ), ANP([QQ(9,1)], h, QQ)]

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(4,1)], h, QQ), [
            ([ANP([QQ(1,1)], h, QQ), ANP([-QQ(3,2), QQ(0,1)], h, QQ)], 1),
            ([ANP([QQ(1,1)], h, QQ), ANP([ QQ(3,2), QQ(0,1)], h, QQ)], 1),
        ])

    f = [ANP([QQ(4,1)], h, QQ), ANP([QQ(8,1)], h, QQ), ANP([QQ(77,1)], h, QQ), ANP([QQ(18,1)], h, QQ), ANP([QQ(153,1)], h, QQ)]

    assert dup_ext_factor(f, K) == \
        (ANP([QQ(4,1)], h, QQ), [
            ([ANP([QQ(1,1)], h, QQ), ANP([-QQ(4,1), QQ(1,1)], h, QQ)], 1),
            ([ANP([QQ(1,1)], h, QQ), ANP([-QQ(3,2), QQ(0,1)], h, QQ)], 1),
            ([ANP([QQ(1,1)], h, QQ), ANP([ QQ(3,2), QQ(0,1)], h, QQ)], 1),
            ([ANP([QQ(1,1)], h, QQ), ANP([ QQ(4,1), QQ(1,1)], h, QQ)], 1),
        ])

def test_dmp_ext_factor():
    h = [QQ(1),QQ(0),QQ(-2)]
    K = QQ.algebraic_field(sqrt(2))

    assert dmp_ext_factor([], 0, K) == (ANP([], h, QQ), [])
    assert dmp_ext_factor([[]], 1, K) == (ANP([], h, QQ), [])

    f = [[ANP([QQ(1)], h, QQ)], [ANP([QQ(1)], h, QQ)]]

    assert dmp_ext_factor(f, 1, K) == (ANP([QQ(1)], h, QQ), [(f, 1)])

    g = [[ANP([QQ(2)], h, QQ)], [ANP([QQ(2)], h, QQ)]]

    assert dmp_ext_factor(g, 1, K) == (ANP([QQ(2)], h, QQ), [(f, 1)])

    f = [[ANP([QQ(1)], h, QQ)], [], [ANP([QQ(-2)], h, QQ), ANP([], h, QQ), ANP([], h, QQ)]]

    assert dmp_ext_factor(f, 1, K) == \
        (ANP([QQ(1)], h, QQ), [
            ([[ANP([QQ(1)], h, QQ)], [ANP([QQ(-1),QQ(0)], h, QQ), ANP([], h, QQ)]], 1),
            ([[ANP([QQ(1)], h, QQ)], [ANP([QQ( 1),QQ(0)], h, QQ), ANP([], h, QQ)]], 1),
        ])

    f = [[ANP([QQ(2)], h, QQ)], [], [ANP([QQ(-4)], h, QQ), ANP([], h, QQ), ANP([], h, QQ)]]

    assert dmp_ext_factor(f, 1, K) == \
        (ANP([QQ(2)], h, QQ), [
            ([[ANP([QQ(1)], h, QQ)], [ANP([QQ(-1),QQ(0)], h, QQ), ANP([], h, QQ)]], 1),
            ([[ANP([QQ(1)], h, QQ)], [ANP([QQ( 1),QQ(0)], h, QQ), ANP([], h, QQ)]], 1),
        ])

def test_dup_factor_list():
    assert dup_factor_list([], ZZ) == (ZZ(0), [])
    assert dup_factor_list([], QQ) == (QQ(0), [])
    assert dup_factor_list([], ZZ['y']) == (DMP([],ZZ), [])
    assert dup_factor_list([], QQ['y']) == (DMP([],QQ), [])

    assert dup_factor_list_include([], ZZ) == [([], 1)]

    assert dup_factor_list([ZZ(7)], ZZ) == (ZZ(7), [])
    assert dup_factor_list([QQ(1,7)], QQ) == (QQ(1,7), [])
    assert dup_factor_list([DMP([ZZ(7)],ZZ)], ZZ['y']) == (DMP([ZZ(7)],ZZ), [])
    assert dup_factor_list([DMP([QQ(1,7)],QQ)], QQ['y']) == (DMP([QQ(1,7)],QQ), [])

    assert dup_factor_list_include([ZZ(7)], ZZ) == [([ZZ(7)], 1)]

    assert dup_factor_list([ZZ(1),ZZ(2),ZZ(1)], ZZ) == \
        (ZZ(1), [([ZZ(1), ZZ(1)], 2)])
    assert dup_factor_list([QQ(1,2),QQ(1),QQ(1,2)], QQ) == \
        (QQ(1,2), [([QQ(1),QQ(1)], 2)])

    assert dup_factor_list_include([ZZ(1),ZZ(2),ZZ(1)], ZZ) == \
        [([ZZ(1), ZZ(1)], 2)]

    K = FF(2)

    assert dup_factor_list([K(1),K(0),K(1)], K) == \
        (K(1), [([K(1), K(1)], 2)])

    assert dup_factor_list([RR(1.0),RR(2.0),RR(1.0)], RR) == \
        (RR(1.0), [([RR(1.0),RR(1.0)], 2)])
    assert dup_factor_list([RR(2.0),RR(4.0),RR(2.0)], RR) == \
        (RR(2.0), [([RR(1.0),RR(1.0)], 2)])

    f = [DMP([ZZ(4),ZZ(0)],ZZ),DMP([ZZ(4),ZZ(0),ZZ(0)],ZZ),DMP([],ZZ)]

    assert dup_factor_list(f, ZZ['y']) == \
        (DMP([ZZ(4)],ZZ), [([DMP([ZZ(1)],ZZ),DMP([],ZZ)], 1),
                           ([DMP([ZZ(1),ZZ(0)],ZZ)], 1),
                           ([DMP([ZZ(1)],ZZ),DMP([ZZ(1),ZZ(0)],ZZ)], 1)])

    f = [DMP([QQ(1,2),QQ(0)],ZZ),DMP([QQ(1,2),QQ(0),QQ(0)],ZZ),DMP([],ZZ)]

    assert dup_factor_list(f, QQ['y']) == \
        (DMP([QQ(1,2)],QQ), [([DMP([QQ(1)],QQ),DMP([],QQ)], 1),
                             ([DMP([QQ(1),QQ(0)],QQ)], 1),
                             ([DMP([QQ(1)],QQ),DMP([QQ(1),QQ(0)],QQ)], 1)])

    raises(DomainError, "dup_factor_list([EX(sin(1))], EX)")

def test_dmp_factor_list():
    assert dmp_factor_list([[]], 1, ZZ) == (ZZ(0), [])
    assert dmp_factor_list([[]], 1, QQ) == (QQ(0), [])
    assert dmp_factor_list([[]], 1, ZZ['y']) == (DMP([],ZZ), [])
    assert dmp_factor_list([[]], 1, QQ['y']) == (DMP([],QQ), [])

    assert dmp_factor_list_include([[]], 1, ZZ) == [([[]], 1)]

    assert dmp_factor_list([[ZZ(7)]], 1, ZZ) == (ZZ(7), [])
    assert dmp_factor_list([[QQ(1,7)]], 1, QQ) == (QQ(1,7), [])
    assert dmp_factor_list([[DMP([ZZ(7)],ZZ)]], 1, ZZ['y']) == (DMP([ZZ(7)],ZZ), [])
    assert dmp_factor_list([[DMP([QQ(1,7)],QQ)]], 1, QQ['y']) == (DMP([QQ(1,7)],QQ), [])

    assert dmp_factor_list_include([[ZZ(7)]], 1, ZZ) == [([[ZZ(7)]], 1)]

    f, g = [ZZ(1),ZZ(2),ZZ(1)], [ZZ(1),ZZ(1)]

    assert dmp_factor_list(dmp_nest(f, 200, ZZ), 200, ZZ) == \
        (ZZ(1), [(dmp_nest(g, 200, ZZ), 2)])

    assert dmp_factor_list(dmp_raise(f, 200, 0, ZZ), 200, ZZ) == \
        (ZZ(1), [(dmp_raise(g, 200, 0, ZZ), 2)])

    assert dmp_factor_list([ZZ(1),ZZ(2),ZZ(1)], 0, ZZ) == \
        (ZZ(1), [([ZZ(1), ZZ(1)], 2)])
    assert dmp_factor_list([QQ(1,2),QQ(1),QQ(1,2)], 0, QQ) == \
        (QQ(1,2), [([QQ(1),QQ(1)], 2)])

    assert dmp_factor_list([[ZZ(1)],[ZZ(2)],[ZZ(1)]], 1, ZZ) == \
        (ZZ(1), [([[ZZ(1)], [ZZ(1)]], 2)])
    assert dmp_factor_list([[QQ(1,2)],[QQ(1)],[QQ(1,2)]], 1, QQ) == \
        (QQ(1,2), [([[QQ(1)],[QQ(1)]], 2)])

    f = [[ZZ(4),ZZ(0)],[ZZ(4),ZZ(0),ZZ(0)],[]]

    assert dmp_factor_list(f, 1, ZZ) == \
        (ZZ(4), [([[ZZ(1)],[]], 1),
                 ([[ZZ(1),ZZ(0)]], 1),
                 ([[ZZ(1)],[ZZ(1),ZZ(0)]], 1)])

    assert dmp_factor_list_include(f, 1, ZZ) == \
        [([[ZZ(4)],[]], 1),
         ([[ZZ(1),ZZ(0)]], 1),
         ([[ZZ(1)],[ZZ(1),ZZ(0)]], 1)]

    f = [[QQ(1,2),QQ(0)],[QQ(1,2),QQ(0),QQ(0)],[]]

    assert dmp_factor_list(f, 1, QQ) == \
        (QQ(1,2), [([[QQ(1)],[]], 1),
                   ([[QQ(1),QQ(0)]], 1),
                   ([[QQ(1)],[QQ(1),QQ(0)]], 1)])

    f = [[RR(2.0)],[],[-RR(8.0),RR(0.0),RR(0.0)]]

    assert dmp_factor_list(f, 1, RR) == \
        (RR(2.0), [([[RR(1.0)],[-RR(2.0),RR(0.0)]], 1),
                   ([[RR(1.0)],[ RR(2.0),RR(0.0)]], 1)])

    f = [[DMP([ZZ(4),ZZ(0)],ZZ)],[DMP([ZZ(4),ZZ(0),ZZ(0)],ZZ)],[DMP([],ZZ)]]

    assert dmp_factor_list(f, 1, ZZ['y']) == \
        (DMP([ZZ(4)],ZZ), [([[DMP([ZZ(1)],ZZ)],[]], 1),
                           ([[DMP([ZZ(1),ZZ(0)],ZZ)]], 1),
                           ([[DMP([ZZ(1)],ZZ)],[DMP([ZZ(1),ZZ(0)],ZZ)]], 1)])

    f = [[DMP([QQ(1,2),QQ(0)],ZZ)],[DMP([QQ(1,2),QQ(0),QQ(0)],ZZ)],[DMP([],ZZ)]]

    assert dmp_factor_list(f, 1, QQ['y']) == \
        (DMP([QQ(1,2)],QQ), [([[DMP([QQ(1)],QQ)],[]], 1),
                             ([[DMP([QQ(1),QQ(0)],QQ)]], 1),
                             ([[DMP([QQ(1)],QQ)],[DMP([QQ(1),QQ(0)],QQ)]], 1)])

    K = FF(2)

    raises(DomainError, "dmp_factor_list([[K(1)],[],[K(1),K(0),K(0)]], 1, K)")
    raises(DomainError, "dmp_factor_list([[EX(sin(1))]], 1, EX)")

