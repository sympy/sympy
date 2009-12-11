from sympy.polys.galoispolys import (
    gf_degree, gf_strip, gf_normal,
    gf_from_dict, gf_to_dict,
    gf_from_int_poly, gf_to_int_poly,
    gf_neg, gf_add_const, gf_sub_const, gf_mul_const, gf_div_const,
    gf_add, gf_sub, gf_add_mul, gf_sub_mul, gf_mul, gf_sqr,
    gf_div, gf_quo, gf_rem, gf_lshift, gf_rshift,
    gf_pow, gf_pow_mod,
    gf_gcd, gf_gcdex,
    gf_LC, gf_TC, gf_monic,
    gf_eval, gf_multi_eval,
    gf_compose, gf_compose_mod,
    gf_trace_map,
    gf_diff, gf_random,
    gf_sqf, gf_sqf_p,
    gf_ddf_zassenhaus, gf_edf_zassenhaus,
    gf_ddf_shoup, gf_edf_shoup,
    gf_factor, gf_factor_sqf)

from sympy import pi, nextprime, raises

def test_gf_degree():
    assert gf_degree([]) == -1
    assert gf_degree([1]) == 0
    assert gf_degree([1,0]) == 1
    assert gf_degree([1,0,0,0,1]) == 4

def test_gf_strip():
    assert gf_strip([]) == []
    assert gf_strip([0]) == []
    assert gf_strip([0,0,0]) == []

    assert gf_strip([1]) == [1]
    assert gf_strip([0,1]) == [1]
    assert gf_strip([0,0,0,1]) == [1]

    assert gf_strip([1,2,0]) == [1,2,0]
    assert gf_strip([0,1,2,0]) == [1,2,0]
    assert gf_strip([0,0,0,1,2,0]) == [1,2,0]

def test_gf_normal():
    assert gf_normal([], 11) == []
    assert gf_normal([1], 11) == [1]
    assert gf_normal([22], 11) == []
    assert gf_normal([12], 11) == [1]

    assert gf_normal([11,22,17,1,0], 11) == [6,1,0]
    assert gf_normal([12,23,17,1,0], 11) == [1,1,6,1,0]

def test_gf_from_to_dict():
    f = {11: 12, 6: 2, 0: 25}
    F = {11: 1, 6: 2, 0: 3}
    g = [1,0,0,0,0,2,0,0,0,0,0,3]

    assert gf_from_dict(f, 11) == g
    assert gf_to_dict(g, 11) == F

    f = {11: -5, 4: 0, 3: 1, 0: 12}
    F = {11: -5, 3: 1, 0: 1}
    g = [6,0,0,0,0,0,0,0,1,0,0,1]

    assert gf_from_dict(f, 11) == g
    assert gf_to_dict(g, 11) == F

def test_gf_from_to_int_poly():
    assert gf_from_int_poly([1,0,7,2,20], 5) == [1,0,2,2,0]
    assert gf_to_int_poly([1,0,4,2,3], 5) == [1,0,-1,2,-2]

def test_gf_monic():
    assert gf_LC([]) == 0
    assert gf_TC([]) == 0

    assert gf_LC([1]) == 1
    assert gf_TC([1]) == 1

    assert gf_LC([1,2]) == 1
    assert gf_TC([1,2]) == 2

    assert gf_monic([], 11) == (0, [])

    assert gf_monic([1], 11) == (1, [1])
    assert gf_monic([2], 11) == (2, [1])

    assert gf_monic([1,2,3,4], 11) == (1, [1,2,3,4])
    assert gf_monic([2,3,4,5], 11) == (2, [1,7,2,8])

def test_gf_arith():
    assert gf_neg([], 11) == []
    assert gf_neg([1], 11) == [10]
    assert gf_neg([1,2,3], 11) == [10,9,8]

    assert gf_add_const([], 0, 11) == []
    assert gf_sub_const([], 0, 11) == []

    assert gf_add_const([], 3, 11) == [3]
    assert gf_sub_const([], 3, 11) == [8]

    assert gf_add_const([1], 3, 11) == [4]
    assert gf_sub_const([1], 3, 11) == [9]

    assert gf_add_const([8], 3, 11) == []
    assert gf_sub_const([3], 3, 11) == []

    assert gf_add_const([1,2,3], 3, 11) == [1,2,6]
    assert gf_sub_const([1,2,3], 3, 11) == [1,2,0]

    assert gf_mul_const([], 0, 11) == []
    assert gf_mul_const([], 1, 11) == []

    assert gf_mul_const([1], 0, 11) == []
    assert gf_mul_const([1], 1, 11) == [1]

    assert gf_mul_const([1,2,3], 0, 11) == []
    assert gf_mul_const([1,2,3], 1, 11) == [1,2,3]
    assert gf_mul_const([1,2,3], 7, 11) == [7,3,10]

    assert gf_add([], [], 11) == []
    assert gf_add([1], [], 11) == [1]
    assert gf_add([], [1], 11) == [1]
    assert gf_add([1], [1], 11) == [2]
    assert gf_add([1], [2], 11) == [3]

    assert gf_add([1,2], [1], 11) == [1,3]
    assert gf_add([1], [1,2], 11) == [1,3]

    assert gf_add([1,2,3], [8,9,10], 11) == [9,0,2]

    assert gf_sub([], [], 11) == []
    assert gf_sub([1], [], 11) == [1]
    assert gf_sub([], [1], 11) == [10]
    assert gf_sub([1], [1], 11) == []
    assert gf_sub([1], [2], 11) == [10]

    assert gf_sub([1,2], [1], 11) == [1,1]
    assert gf_sub([1], [1,2], 11) == [10,10]

    assert gf_sub([3,2,1], [8,9,10], 11) == [6,4,2]

    assert gf_add_mul([1,5,6], [7,3], [8,0,6,1], 11) == [1,2,10,8,9]
    assert gf_sub_mul([1,5,6], [7,3], [8,0,6,1], 11) == [10,9,3,2,3]

    assert gf_mul([], [], 11) == []
    assert gf_mul([], [1], 11) == []
    assert gf_mul([1], [], 11) == []
    assert gf_mul([1], [1], 11) == [1]
    assert gf_mul([5], [7], 11) == [2]

    assert gf_mul([3,0,0,6,1,2], [4,0,1,0], 11) == [1,0,3,2,4,3,1,2,0]
    assert gf_mul([4,0,1,0], [3,0,0,6,1,2], 11) == [1,0,3,2,4,3,1,2,0]

    assert gf_mul([2,0,0,1,7], [2,0,0,1,7], 11) == [4,0,0,4,6,0,1,3,5]

    assert gf_sqr([], 11) == []
    assert gf_sqr([2], 11) == [4]
    assert gf_sqr([1,2], 11) == [1,4,4]

    assert gf_sqr([2,0,0,1,7], 11) == [4,0,0,4,6,0,1,3,5]

def test_gf_division():
    raises(ZeroDivisionError, "gf_div([1,2,3], [], 11)")
    raises(ZeroDivisionError, "gf_quo([1,2,3], [], 11)")
    raises(ZeroDivisionError, "gf_rem([1,2,3], [], 11)")

    f, g, q, r = [5,4,3,2,1,0], [1,2,3], [5,1,0,6], [3,3]

    assert gf_div(f, g, 7) == (q, r)
    assert gf_quo(f, g, 7) == q
    assert gf_rem(f, g, 7) == r

    f, g, q, r = [5,4,3,2,1,0], [1,2,3,0], [5,1,0], [6,1,0]

    assert gf_div(f, g, 7) == (q, r)
    assert gf_quo(f, g, 7) == q
    assert gf_rem(f, g, 7) == r

def test_gf_shift():
    f = [1,2,3,4,5]

    assert gf_lshift([], 5) == []
    assert gf_rshift([], 5) == ([], [])

    assert gf_lshift(f, 1) == [1,2,3,4,5,0]
    assert gf_lshift(f, 2) == [1,2,3,4,5,0,0]

    assert gf_rshift(f, 0) == (f, [])
    assert gf_rshift(f, 1) == ([1,2,3,4], [5])
    assert gf_rshift(f, 3) == ([1,2], [3,4,5])
    assert gf_rshift(f, 5) == ([], f)

def test_gf_powering():
    assert gf_pow([1,0,0,1,8], 0, 11) == [1]
    assert gf_pow([1,0,0,1,8], 1, 11) == [1, 0, 0, 1, 8]
    assert gf_pow([1,0,0,1,8], 2, 11) == [1, 0, 0, 2, 5, 0, 1, 5, 9]

    assert gf_pow([1,0,0,1,8], 5, 11) == \
        [1, 0, 0, 5, 7, 0, 10, 6, 2, 10, 9, 6, 10, 6, 6, 0, 5, 2, 5, 9, 10]

    assert gf_pow([1,0,0,1,8], 8, 11) == \
        [1, 0, 0, 8, 9, 0, 6, 8, 10, 1, 2, 5, 10, 7, 7, 9, 1, 2, 0, 0, 6, 2,
         5, 2, 5, 7, 7, 9, 10, 10, 7, 5, 5]

    assert gf_pow([1,0,0,1,8], 45, 11) == \
        [ 1, 0, 0,  1, 8, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0,
          0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  4, 0, 0,  4, 10, 0, 0, 0, 0, 0, 0,
         10, 0, 0, 10, 3, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0,  0, 0, 0, 0, 0, 0, 0,
          6, 0, 0,  6, 4, 0, 0, 0, 0, 0, 0,  8, 0, 0,  8,  9, 0, 0, 0, 0, 0, 0,
         10, 0, 0, 10, 3, 0, 0, 0, 0, 0, 0,  4, 0, 0,  4, 10, 0, 0, 0, 0, 0, 0,
          8, 0, 0,  8, 9, 0, 0, 0, 0, 0, 0,  9, 0, 0,  9,  6, 0, 0, 0, 0, 0, 0,
          3, 0, 0,  3, 2, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10,  3, 0, 0, 0, 0, 0, 0,
         10, 0, 0, 10, 3, 0, 0, 0, 0, 0, 0,  2, 0, 0,  2,  5, 0, 0, 0, 0, 0, 0,
          4, 0, 0, 4, 10]

    assert gf_pow_mod([1,0,0,1,8], 0, [2,0,7], 11) == [1]
    assert gf_pow_mod([1,0,0,1,8], 1, [2,0,7], 11) == [1,1]
    assert gf_pow_mod([1,0,0,1,8], 2, [2,0,7], 11) == [2,3]
    assert gf_pow_mod([1,0,0,1,8], 5, [2,0,7], 11) == [7,8]
    assert gf_pow_mod([1,0,0,1,8], 8, [2,0,7], 11) == [1,5]
    assert gf_pow_mod([1,0,0,1,8], 45, [2,0,7], 11) == [5,4]

def test_gf_euclidean():
    assert gf_gcd([], [], 11) == []
    assert gf_gcd([2], [], 11) == [1]
    assert gf_gcd([], [2], 11) == [1]
    assert gf_gcd([2], [2], 11) == [1]

    assert gf_gcd([], [1,0], 11) == [1,0]
    assert gf_gcd([1,0], [], 11) == [1,0]

    assert gf_gcd([3,0], [3,0], 11) == [1,0]

    assert gf_gcd([1,8,7], [1,7,1,7], 11) == [1,7]

    assert gf_gcdex([], [], 11) == ([1], [], [])
    assert gf_gcdex([2], [], 11) == ([6], [], [1])
    assert gf_gcdex([], [2], 11) == ([], [6], [1])
    assert gf_gcdex([2], [2], 11) == ([], [6], [1])

    assert gf_gcdex([], [3,0], 11) == ([], [4], [1,0])
    assert gf_gcdex([3,0], [], 11) == ([4], [], [1,0])

    assert gf_gcdex([3,0], [3,0], 11) == ([], [4], [1,0])

    assert gf_gcdex([1,8,7], [1,7,1,7], 11) == ([5,6], [6], [1,7])

def test_gf_diff():
    assert gf_diff([], 11) == []
    assert gf_diff([7], 11) == []

    assert gf_diff([7,3], 11) == [7]
    assert gf_diff([7,3,1], 11) == [3,3]

    assert gf_diff([1,0,0,0,0,0,0,0,0,0,0,1], 11) == []

def test_gf_eval():
    assert gf_eval([], 4, 11) == 0
    assert gf_eval([], 27, 11) == 0
    assert gf_eval([7], 4, 11) == 7
    assert gf_eval([7], 27, 11) == 7

    assert gf_eval([1,0,3,2,4,3,1,2,0], 0, 11) == 0
    assert gf_eval([1,0,3,2,4,3,1,2,0], 4, 11) == 9
    assert gf_eval([1,0,3,2,4,3,1,2,0], 27, 11) == 5

    assert gf_eval([4,0,0,4,6,0,1,3,5], 0, 11) == 5
    assert gf_eval([4,0,0,4,6,0,1,3,5], 4, 11) == 3
    assert gf_eval([4,0,0,4,6,0,1,3,5], 27, 11) == 9

    assert gf_multi_eval([3,2,1], [0,1,2,3], 11) == [1,6,6,1]

def test_gf_compose():
    assert gf_compose([], [1,0], 11) == []

    assert gf_compose([1], [], 11) == [1]
    assert gf_compose([1,0], [], 11) == []
    assert gf_compose([1,0], [1,0], 11) == [1,0]

    f, g, h = [1, 1, 4, 9, 1], [1,1,1], [1,0,0,2]

    assert gf_compose(g, h, 11) == [1,0,0,5,0,0,7]
    assert gf_compose_mod(g, h, f, 11) == [3,9,6,10]

def test_gf_trace_map():
    f, a, c = [1, 1, 4, 9, 1], [1,1,1], [1,0]
    b = gf_pow_mod(c, 11, f, 11)

    assert gf_trace_map(a, b, c, 0, f, 11) == \
        ([1, 1, 1], [1, 1, 1])
    assert gf_trace_map(a, b, c, 1, f, 11) == \
        ([5, 2, 10, 3], [5, 3, 0, 4])
    assert gf_trace_map(a, b, c, 2, f, 11) == \
        ([5, 9, 5, 3], [10, 1, 5, 7])
    assert gf_trace_map(a, b, c, 3, f, 11) == \
        ([1, 10, 6, 0], [7])
    assert gf_trace_map(a, b, c, 4, f, 11) == \
        ([1, 1, 1], [1, 1, 8])
    assert gf_trace_map(a, b, c, 5, f, 11) == \
        ([5, 2, 10, 3], [5, 3, 0, 0])
    assert gf_trace_map(a, b, c, 11, f, 11) == \
        ([1, 10, 6, 0], [10])

def test_gf_squarefree():
    assert gf_sqf([], 11) == (0, [])
    assert gf_sqf([1], 11) == (1, [])
    assert gf_sqf([1,1], 11) == (1, [([1, 1], 1)])

    assert gf_sqf_p([], 11) == True
    assert gf_sqf_p([1], 11) == True
    assert gf_sqf_p([1,1], 11) == True

    f = gf_from_dict({11: 1, 0: 1}, 11)

    assert gf_sqf_p(f, 11) == False

    assert gf_sqf(f, 11) == \
       (1, [([1, 1], 11)])

    f = [1, 5, 8, 4]

    assert gf_sqf_p(f, 11) == False

    assert gf_sqf(f, 11) == \
        (1, [([1, 1], 1),
             ([1, 2], 2)])

def test_gf_ddf():
    f = gf_from_dict({15: 1, 0: -1}, 11)
    g = [([1, 0, 0, 0, 0, 10], 1),
         ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], 2)]

    assert gf_ddf_zassenhaus(f, 11) == g
    assert gf_ddf_shoup(f, 11) == g

    f = gf_from_dict({63: 1, 0: 1}, 2)
    g = [([1, 1], 1),
         ([1, 1, 1], 2),
         ([1, 1, 1, 1, 1, 1, 1], 3),
         ([1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1], 6)]

    assert gf_ddf_zassenhaus(f, 2) == g
    assert gf_ddf_shoup(f, 2) == g

    f = gf_from_dict({6: 1, 5: -1, 4: 1, 3: 1, 1: -1}, 3)
    g = [([1, 1, 0], 1),
         ([1, 1, 0, 1, 2], 2)]

    assert gf_ddf_zassenhaus(f, 3) == g
    assert gf_ddf_shoup(f, 3) == g

    f = [1, 2, 5, 26, 677, 436, 791, 325, 456, 24, 577]
    g = [([1, 701], 1),
         ([1, 110, 559, 532, 694, 151, 110, 70, 735, 122], 9)]

    assert gf_ddf_zassenhaus(f, 809) == g
    assert gf_ddf_shoup(f, 809) == g

    p = nextprime(int((2**15 * pi).evalf()))
    f = gf_from_dict({15: 1, 1: 1, 0: 1}, p)
    g = [([1, 22730, 68144], 2),
         ([1, 64876, 83977, 10787, 12561, 68608, 52650, 88001, 84356], 4),
         ([1, 15347, 95022, 84569, 94508, 92335], 5)]

    assert gf_ddf_zassenhaus(f, p) == g
    assert gf_ddf_shoup(f, p) == g

def test_gf_edf():
    f = [1, 1, 0, 1, 2]
    g = [[1, 0, 1], [1, 1, 2]]

    assert gf_edf_zassenhaus(f, 2, 3) == g
    assert gf_edf_shoup(f, 2, 3) == g

def test_gf_factor():
    assert gf_factor([], 11) == (0, [])
    assert gf_factor([1], 11) == (1, [])
    assert gf_factor([1,1], 11) == (1, [([1, 1], 1)])

    f, p = [1,0,0,1,0], 2

    g = (1, [([1, 0], 1),
             ([1, 1], 1),
             ([1, 1, 1], 1)])

    assert gf_factor(f, p, method='zassenhaus') == g
    assert gf_factor(f, p, method='shoup') == g

    g = (1, [[1, 0],
             [1, 1],
             [1, 1, 1]])

    assert gf_factor_sqf(f, p, method='zassenhaus') == g
    assert gf_factor_sqf(f, p, method='shoup') == g

    assert gf_factor([1, 5, 8, 4], 11) == \
        (1, [([1, 1], 1), ([1, 2], 2)])

    assert gf_factor([1, 1, 10, 1, 0, 10, 10, 10, 0, 0], 11) == \
        (1, [([1, 0], 2), ([1, 9, 5], 1), ([1, 3, 0, 8, 5, 2], 1)])

    assert gf_factor(gf_from_dict({32: 1, 0: 1}, 11), 11) == \
        (1, [([1, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 10], 1),
             ([1, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 10], 1)])

    assert gf_factor(gf_from_dict({32: 8, 0: 5}, 11), 11) == \
        (8, [([1, 3], 1),
             ([1, 8], 1),
             ([1, 0, 9], 1),
             ([1, 2, 2], 1),
             ([1, 9, 2], 1),
             ([1, 0, 5, 0, 7], 1),
             ([1, 0, 6, 0, 7], 1),
             ([1, 0, 0, 0, 1, 0, 0, 0, 6], 1),
             ([1, 0, 0, 0, 10, 0, 0, 0, 6], 1)])

    assert gf_factor(gf_from_dict({63: 8, 0: 5}, 11), 11) == \
        (8, [([1, 7], 1),
             ([1, 4, 5], 1),
             ([1, 6, 8, 2], 1),
             ([1, 9, 9, 2], 1),
             ([1, 0, 0, 9, 0, 0, 4], 1),
             ([1, 2, 0, 8, 4, 6, 4], 1),
             ([1, 2, 3, 8, 0, 6, 4], 1),
             ([1, 2, 6, 0, 8, 4, 4], 1),
             ([1, 3, 3, 1, 6, 8, 4], 1),
             ([1, 5, 6, 0, 8, 6, 4], 1),
             ([1, 6, 2, 7, 9, 8, 4], 1),
             ([1, 10, 4, 7, 10, 7, 4], 1),
             ([1, 10, 10, 1, 4, 9, 4], 1)])

    # Gathen polynomials: x**n + x + 1 (mod p > 2**n * pi)

    p = nextprime(int((2**15 * pi).evalf()))
    f = gf_from_dict({15: 1, 1: 1, 0: 1}, p)

    assert gf_sqf_p(f, p) == True

    g = (1, [([1, 22730, 68144], 1),
             ([1, 81553, 77449, 86810, 4724], 1),
             ([1, 86276, 56779, 14859, 31575], 1),
             ([1, 15347, 95022, 84569, 94508, 92335], 1)])

    assert gf_factor(f, p, method='zassenhaus') == g
    assert gf_factor(f, p, method='shoup') == g

    g = (1, [[1, 22730, 68144],
             [1, 81553, 77449, 86810, 4724],
             [1, 86276, 56779, 14859, 31575],
             [1, 15347, 95022, 84569, 94508, 92335]])

    assert gf_factor_sqf(f, p, method='zassenhaus') == g
    assert gf_factor_sqf(f, p, method='shoup') == g

    # Shoup polynomials: f = a_0 x**n + a_1 x**(n-1) + ... + a_n
    # (mod p > 2**(n-2) * pi), where a_n = a_{n-1}**2 + 1, a_0 = 1

    p = nextprime(int((2**4 * pi).evalf()))
    f = [1, 2, 5, 26, 41, 39, 38]  # deg(f) = 6

    assert gf_sqf_p(f, p) == True

    g = (1, [([1, 44, 26], 1),
             ([1, 11, 25, 18, 30], 1)])

    assert gf_factor(f, p, method='zassenhaus') == g
    assert gf_factor(f, p, method='shoup') == g

    g = (1, [[1, 44, 26],
             [1, 11, 25, 18, 30]])

    assert gf_factor_sqf(f, p, method='zassenhaus') == g
    assert gf_factor_sqf(f, p, method='shoup') == g
