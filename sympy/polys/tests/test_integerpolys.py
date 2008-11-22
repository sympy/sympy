
from sympy.polys.integerpolys import (
    zzx_degree, zzx_strip,
    zzx_from_dict, zzx_to_dict,
    zzx_LC, zzx_TC,
    zzx_abs, zzx_neg,
    zzx_add_term, zzx_sub_term, zzx_mul_term,
    zzx_add, zzx_sub, zzx_mul, zzx_sqr,
    zzx_div, zzx_quo, zzx_rem,
    zzx_heu_gcd, zzx_mod_gcd,
    zzx_max_norm, zzx_l1_norm,
    zzx_diff, zzx_eval, zzx_trunc,
    zzx_sqf_part,
    zzx_content, zzx_primitive,
    zzx_hensel_step, zzx_hensel_lift,
    zzx_zassenhaus, zzx_factor)

from sympy import raises

def test_zzx_degree():
    assert zzx_degree([]) == -1
    assert zzx_degree([1]) == 0
    assert zzx_degree([1,0]) == 1
    assert zzx_degree([1,0,0,0,1]) == 4

def test_zzx_strip():
    assert zzx_strip([]) == []
    assert zzx_strip([0]) == []
    assert zzx_strip([0,0,0]) == []

    assert zzx_strip([1]) == [1]
    assert zzx_strip([0,1]) == [1]
    assert zzx_strip([0,0,0,1]) == [1]

    assert zzx_strip([1,2,0]) == [1,2,0]
    assert zzx_strip([0,1,2,0]) == [1,2,0]
    assert zzx_strip([0,0,0,1,2,0]) == [1,2,0]

def test_zzx_from_to_dict():
    f = {8: 3, 5: 2, 0: 8}
    g = [3,0,0,2,0,0,0,0,8]

    assert zzx_from_dict(f) == g
    assert zzx_to_dict(g) == f

def test_zzx_add_term():
    assert zzx_add_term([], 0, 0) == []

    assert zzx_add_term([], 1, 0) == [1]
    assert zzx_add_term([], 1, 1) == [1, 0]
    assert zzx_add_term([], 1, 2) == [1, 0, 0]

    assert zzx_add_term([1,1,1], 1, 0) == [1, 1, 2]
    assert zzx_add_term([1,1,1], 1, 1) == [1, 2, 1]
    assert zzx_add_term([1,1,1], 1, 2) == [2, 1, 1]

    assert zzx_add_term([1,1,1], 1, 3) == [1, 1, 1, 1]
    assert zzx_add_term([1,1,1], 1, 4) == [1, 0, 1, 1, 1]
    assert zzx_add_term([1,1,1], 1, 5) == [1, 0, 0, 1, 1, 1]
    assert zzx_add_term([1,1,1], 1, 6) == [1, 0, 0, 0, 1, 1, 1]

    assert zzx_add_term([1,1,1], -1, 2) == [1, 1]

def test_zzx_sub_term():
    assert zzx_sub_term([], 0, 0) == []

    assert zzx_sub_term([], 1, 0) == [-1]
    assert zzx_sub_term([], 1, 1) == [-1, 0]
    assert zzx_sub_term([], 1, 2) == [-1, 0, 0]

    assert zzx_sub_term([1,1,1], 2, 0) == [ 1, 1,-1]
    assert zzx_sub_term([1,1,1], 2, 1) == [ 1,-1, 1]
    assert zzx_sub_term([1,1,1], 2, 2) == [-1, 1, 1]

    assert zzx_sub_term([1,1,1], 1, 3) == [-1, 1, 1, 1]
    assert zzx_sub_term([1,1,1], 1, 4) == [-1, 0, 1, 1, 1]
    assert zzx_sub_term([1,1,1], 1, 5) == [-1, 0, 0, 1, 1, 1]
    assert zzx_sub_term([1,1,1], 1, 6) == [-1, 0, 0, 0, 1, 1, 1]

    assert zzx_sub_term([1,1,1], 1, 2) == [1, 1]

def test_zzx_mul_term():
    assert zzx_mul_term([],    2, 3) == []
    assert zzx_mul_term([1,1], 0, 3) == []

    assert zzx_mul_term([1,2,3], 2, 0) == [2,4,6]
    assert zzx_mul_term([1,2,3], 2, 1) == [2,4,6,0]
    assert zzx_mul_term([1,2,3], 2, 2) == [2,4,6,0,0]
    assert zzx_mul_term([1,2,3], 2, 3) == [2,4,6,0,0,0]

def test_zzx_abs():
    assert zzx_abs([-1,2,3]) == [1,2,3]

def test_zzx_neg():
    assert zzx_neg([-1,2,3]) == [1,-2,-3]

def test_zzx_add():
    assert zzx_add([], []) == []
    assert zzx_add([1], []) == [1]
    assert zzx_add([], [1]) == [1]
    assert zzx_add([1], [1]) == [2]
    assert zzx_add([1], [2]) == [3]

    assert zzx_add([1,2], [1]) == [1,3]
    assert zzx_add([1], [1,2]) == [1,3]

    assert zzx_add([1,2,3], [8,9,10]) == [9,11,13]

def test_zzx_sub():
    assert zzx_sub([], []) == []
    assert zzx_sub([1], []) == [1]
    assert zzx_sub([], [1]) == [-1]
    assert zzx_sub([1], [1]) == []
    assert zzx_sub([1], [2]) == [-1]

    assert zzx_sub([1,2], [1]) == [1,1]
    assert zzx_sub([1], [1,2]) == [-1,-1]

    assert zzx_sub([3,2,1], [8,9,10]) == [-5,-7,-9]

def test_zzx_mul():
    assert zzx_mul([], []) == []
    assert zzx_mul([], [1]) == []
    assert zzx_mul([1], []) == []
    assert zzx_mul([1], [1]) == [1]
    assert zzx_mul([5], [7]) == [35]

    assert zzx_mul([3,0,0,6,1,2], [4,0,1,0]) == [12,0,3,24,4,14,1,2,0]
    assert zzx_mul([4,0,1,0], [3,0,0,6,1,2]) == [12,0,3,24,4,14,1,2,0]

    assert zzx_mul([2,0,0,1,7], [2,0,0,1,7]) == [4,0,0,4,28,0,1,14,49]

def test_zzx_sqr():
    assert zzx_sqr([]) == []
    assert zzx_sqr([2]) == [4]
    assert zzx_sqr([1,2]) == [1,4,4]

    assert zzx_sqr([2,0,0,1,7]) == [4,0,0,4,28,0,1,14,49]

def test_zzx_div():
    raises(ZeroDivisionError, "zzx_div([1,2,3], [])")
    raises(ZeroDivisionError, "zzx_quo([1,2,3], [])")
    raises(ZeroDivisionError, "zzx_rem([1,2,3], [])")

    f, g, q, r = [5,4,3,2,1], [1,2,3], [5,-6,0], [20,1]

    assert zzx_div(f, g) == (q, r)
    assert zzx_quo(f, g) == q
    assert zzx_rem(f, g) == r

    f, g, q, r = [5,4,3,2,1,0], [1,2,0,0,9], [5,-6], [15,2,-44,54]

    assert zzx_div(f, g) == (q, r)
    assert zzx_quo(f, g) == q
    assert zzx_rem(f, g) == r

def test_zzx_gcd():
    assert zzx_heu_gcd([], []) == ([], [], [])
    assert zzx_mod_gcd([], []) == ([], [], [])

    assert zzx_heu_gcd([2], []) == ([2], [1], [])
    assert zzx_mod_gcd([2], []) == ([2], [1], [])

    assert zzx_heu_gcd([2], [2]) == ([2], [1], [1])
    assert zzx_mod_gcd([2], [2]) == ([2], [1], [1])

    assert zzx_heu_gcd([1,2,1], [1]) == ([1], [1, 2, 1], [1])
    assert zzx_mod_gcd([1,2,1], [1]) == ([1], [1, 2, 1], [1])

    assert zzx_heu_gcd([1,2,1], [2]) == ([1], [1, 2, 1], [2])
    assert zzx_mod_gcd([1,2,1], [2]) == ([1], [1, 2, 1], [2])

    assert zzx_heu_gcd([2,4,2], [2]) == ([2], [1, 2, 1], [1])
    assert zzx_mod_gcd([2,4,2], [2]) == ([2], [1, 2, 1], [1])

    assert zzx_heu_gcd([2], [2,4,2]) == ([2], [1], [1, 2, 1])
    assert zzx_mod_gcd([2], [2,4,2]) == ([2], [1], [1, 2, 1])

    assert zzx_heu_gcd([2,4,2], [1,1]) == ([1, 1], [2, 2], [1])
    assert zzx_mod_gcd([2,4,2], [1,1]) == ([1, 1], [2, 2], [1])

    assert zzx_heu_gcd([1,1], [2,4,2]) == ([1, 1], [1], [2, 2])
    assert zzx_mod_gcd([1,1], [2,4,2]) == ([1, 1], [1], [2, 2])

    f = [1,8,21,22,8]
    g = [1,6,11,6]

    h = [1,3,2]

    cff = [1,5,4]
    cfg = [1,3]

    assert zzx_heu_gcd(f, g) == (h, cff, cfg)
    assert zzx_mod_gcd(f, g) == (h, cff, cfg)

    f = [1,0,0,0,-4]
    g = [1,0,4,0, 4]

    h = [1,0,2]

    cff = [1,0,-2]
    cfg = [1,0, 2]

    assert zzx_heu_gcd(f, g) == (h, cff, cfg)
    assert zzx_mod_gcd(f, g) == (h, cff, cfg)

def test_zzx_norm():
    assert zzx_max_norm([]) == 0
    assert zzx_max_norm([1]) == 1
    assert zzx_max_norm([1,4,2,3]) == 4

    assert zzx_l1_norm([]) == 0
    assert zzx_l1_norm([1]) == 1
    assert zzx_l1_norm([1,4,2,3]) == 10

def test_zzx_diff():
    assert zzx_diff([]) == []
    assert zzx_diff([7]) == []
    assert zzx_diff([2,7]) == [2]
    assert zzx_diff([1,2,1]) == [2,2]
    assert zzx_diff([1,2,3,4]) == [3,4,3]
    assert zzx_diff([1,-1,0,0,2]) == [4,-3,0,0]

def test_zzx_eval():
    assert zzx_eval([], 7) == 0
    assert zzx_eval([1,2,3], 7) == 66

def test_zzx_trunc():
    assert zzx_trunc([1,2,3,4,5,6], 3) == [1, -1, 0, 1, -1, 0]
    assert zzx_trunc([6,5,4,3,2,1], 3) ==    [-1, 1, 0, -1, 1]

def test_zzx_content():
    assert zzx_content([]) == 0
    assert zzx_content([1]) == 1
    assert zzx_content([1,1]) == 1
    assert zzx_content([2,2]) == 2
    assert zzx_content([1,2,1]) == 1
    assert zzx_content([2,4,2]) == 2

def test_zzx_primitive():
    assert zzx_primitive([]) == (0, [])
    assert zzx_primitive([1]) == (1, [1])
    assert zzx_primitive([1,1]) == (1, [1,1])
    assert zzx_primitive([2,2]) == (2, [1,1])
    assert zzx_primitive([1,2,1]) == (1, [1,2,1])
    assert zzx_primitive([2,4,2]) == (2, [1,2,1])

def test_zzx_sqf_part():
    assert zzx_sqf_part([2,3,0,0]) == [2,3,0]
    assert zzx_sqf_part([1,0,1,1]) == [1,0,1,1]

def test_zzx_hensel_step():
    f = zzx_from_dict({4:1, 0:-1})

    g = zzx_from_dict({3:1, 2:2, 1:-1, 0:-2})
    h = zzx_from_dict({1:1, 0:-2})
    s = zzx_from_dict({0:-2})
    t = zzx_from_dict({2:2, 1:-2, 0:-1})

    G, H, S, T = zzx_hensel_step(5, f, g, h, s, t)

    assert G == zzx_from_dict({3:1, 2:7, 1:-1, 0:-7})
    assert H == zzx_from_dict({1:1, 0:-7})
    assert S == zzx_from_dict({0:8})
    assert T == zzx_from_dict({2:-8, 1:-12, 0:-1})

def test_zzx_hensel_lift():
    f = zzx_from_dict({4:1, 0:-1})

    f1 = zzx_from_dict({1:1, 0:-1})
    f2 = zzx_from_dict({1:1, 0:-2})
    f3 = zzx_from_dict({1:1, 0: 2})
    f4 = zzx_from_dict({1:1, 0: 1})

    ff_list = zzx_hensel_lift(5, f, [f1, f2, f3, f4], 4)

    assert zzx_to_dict(ff_list[0]) == {0: -1,   1: 1}
    assert zzx_to_dict(ff_list[1]) == {0: -182, 1: 1}
    assert zzx_to_dict(ff_list[2]) == {0:  182, 1: 1}
    assert zzx_to_dict(ff_list[3]) == {0:  1,   1: 1}

def test_zzx_factor():
    f = [1,0,0,1,1]

    for i in xrange(0, 20):
        assert zzx_factor(f) == (1, [(f, 1)])

    assert zzx_factor([2,4,2]) == \
        (2, [([1, 1], 2)])

    assert zzx_factor([1,-6,11,-6]) == \
        (1, [([1,-3], 1),
             ([1,-2], 1),
             ([1,-1], 1)])

    assert zzx_factor([-1,0,0,0,1,0,0]) == \
        (-1, [([1,-1], 1),
              ([1, 1], 1),
              ([1, 0], 2),
              ([1, 0, 1], 1)])

    assert zzx_factor(zzx_from_dict({10:1, 0:-1})) == \
        (1, [([1,-1], 1),
             ([1, 1], 1),
             ([1,-1, 1,-1, 1], 1),
             ([1, 1, 1, 1, 1], 1)])

