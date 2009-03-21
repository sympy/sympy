
from sympy.polys.integerpolys import (
    INT_TYPE,
    poly_LC, poly_TC,
    poly_nth,
    poly_level,
    poly_univariate_p)

from sympy.polys.integerpolys import (
    zzx_degree, zzX_degree,
    zzx_strip, zzX_strip,
    zzX_zz_LC, zzX_zz_TC,
    zzX_zero, zzX_zero_of,
    zzX_const, zzX_const_of,
    zzX_value, zzX_lift,
    zzX_zeros_of,
    zzX_zero_p, zzX_one_p,
    zzx_from_dict, zzX_from_dict,
    zzx_from_poly, zzX_from_poly,
    zzx_to_dict, zzX_to_dict,
    zzx_to_poly, zzX_to_poly,
    zzx_abs, zzX_abs,
    zzx_neg, zzX_neg,
    zzx_add_term, zzX_add_term,
    zzx_sub_term, zzX_sub_term,
    zzx_mul_term, zzX_mul_term,
    zzx_mul_const, zzX_mul_const,
    zzx_quo_const, zzX_quo_const,
    zzx_compose_term,
    zzx_add, zzX_add,
    zzx_sub, zzX_sub,
    zzx_add_mul, zzX_add_mul,
    zzx_sub_mul, zzX_sub_mul,
    zzx_mul, zzX_mul,
    zzx_sqr, zzX_sqr,
    zzx_div, zzX_div,
    zzx_quo, zzX_quo,
    zzx_rem, zzX_rem,
    zzx_max_norm, zzX_max_norm,
    zzx_l1_norm, zzX_l1_norm,
    zzx_diff, zzX_diff,
    zzx_eval, zzX_eval,
    zzX_eval_coeffs,
    zzx_trunc, zzX_trunc,
    zzx_content, zzX_content, zzX_zz_content,
    zzx_primitive, zzX_primitive, zzX_zz_primitive,
    zzx_sqf_part, zzX_sqf_part,
    zzx_sqf_p, zzX_sqf_p,
    zzx_gcd, zzX_gcd,
    zzx_cofactors, zzX_cofactors,
    zzx_heu_gcd, zzX_heu_gcd,
    zzx_mod_gcd,
    zzx_hensel_step, zzx_hensel_lift,
    zzx_zassenhaus, zzx_factor, zzx_cyclotomic_factor)

from sympy.polys.integerpolys import (
    HeuristicGCDFailed, ExactQuotientFailed)

from sympy.polys.specialpolys import (
    zzX_fateman_poly_F_1, zzX_fateman_poly_F_2, zzX_fateman_poly_F_3)

from sympy import raises, symbols

x, y, z, t = symbols('x,y,z,t')

F_0 = (((z**2 + 2*z + 3)*y + 2)*x**2 + 3*x + (4*z**2 + 5*z + 6)*y**2 + (z**2 + 2*z + 1)*y + 1).as_poly(x,y,z)
f_0 = [
    [[1,2,3], [2]],
    [[3]],
    [[4,5,6], [1,2,1], [1]]
]

F_1 = ((z + x*y + 10)*(x*z + y + 30)*(y*z + x + 20)).as_poly(x,y,z)
f_1 = [
    [[1, 0], []],
    [[1, 0, 1], [20, 30], [1, 10, 0]],
    [[1, 0], [30, 20], [1, 10, 1, 610], [20, 230, 300]],
    [[1, 10, 0], [30, 320, 200], [600, 6000]]
]

F_2 = ((x**3*(z + y) + z - 11)*(x**2*(z**2 + y**2) + y + 90)).as_poly(x,y,z)
f_2 = [
    [[1], [1, 0], [1, 0, 0], [1, 0, 0, 0]],
    [[]],
    [[1], [1, 90], [90, 0]],
    [[1, -11], [], [1, -11, 0, 0]],
    [[]],
    [[1, -11], [90, -990]]
]

F_3 = ((y*z**3 + x*y*z + y**2 + x**3)*(x*(z**4 + 1) + z + x**2*y**2)).as_poly(x,y,z)
f_3 = [
    [[1], [], []],
    [[1, 0, 0, 0, 1]],
    [[1, 0], [], [], [1, 0]],
    [[1], [1, 0, 0, 0], [], [1, 0, 0, 0, 1, 0], []],
    [[1, 0, 0, 0, 1], [1, 0, 0, 0, 1, 1, 0, 0], []],
    [[1, 0], [1, 0, 0, 0, 0], []]
]

F_4 = ((z**2 - x**3*y + 3)*(z**2 + x*y**3)*(z**2 + x**3*y**4)*(y**4*z**2 + x**2*z + 5)).as_poly(x,y,z)
f_4 = [
    [[-1, 0], [], [], [], [], [], [], [], []],
    [[-1, 0, 0, 0], [], [], [], [], []],
    [[-1, 0, 0], [], [], [], [-5], [], [], [], [], [], [], [], []],
    [[-1, 0, 0, 0, 0], [], [1, 0, 3, 0], [], [-5, 0, 0], [-1, 0, 0, 0], [], [], [], []],
    [[1, 0, 3, 0, 0, 0], [], [], [-1, 0, 0, 0, 0, 0], []],
    [[1, 0, 3, 0, 0], [], [], [-1, 0, 0, 0, 0], [5, 0, 15], [], [], [-5, 0, 0], [], [], [], []],
    [[1, 0, 3, 0, 0, 0, 0], [], [], [-1, 0, 0, 0, 0, 0, 0], [5, 0, 15, 0, 0], [1, 0, 3, 0, 0, 0], [], [-5, 0, 0, 0, 0], []],
    [[1, 0, 3, 0, 0, 0, 0, 0]],
    [[1, 0, 3, 0, 0, 0, 0], [], [], [], [5, 0, 15, 0, 0], [], [], []],
    [[1, 0, 3, 0, 0, 0, 0, 0, 0], [], [], [], [5, 0, 15, 0, 0, 0, 0]]
]

F_5 = ((z - y - x)**3).as_poly(x,y,z)
f_5 = [
    [[-1]],
    [[-3], [3, 0]],
    [[-3], [6, 0], [-3, 0, 0]],
    [[-1], [3, 0], [-3, 0, 0], [1, 0, 0, 0]]
]

F_6 = ((3*z**3 + 2*t*z - 9*y**3 - y**2 + 45*x**3)*(t**2*z**3 + 47*x*y - t**2)).as_poly(x,y,z,t)
f_6 = [
    [[[2115]], [[]]],
    [[[45, 0, 0], [], [], [-45, 0, 0]]],
    [[[]]],
    [[[-423]], [[-47]], [[]], [[141], [], [94, 0], []], [[]]],
    [[[-9, 0, 0], [], [], [9, 0, 0]],
     [[-1, 0, 0], [], [], [1, 0, 0]],
     [[]],
     [[3, 0, 0], [], [2, 0, 0, 0], [-3, 0, 0], [], [-2, 0, 0, 0], []]
    ]
]

g_0 = [
    [[24], [48], [], []],
    [[24], [], [], [-72], [], []],
    [[25], [2], [], [4], [8]],
    [[1], [], [], [1], [], [], [-12]],
    [[1], [-1], [-2], [292], [], []],
    [[-1], [], [], [3], [], [], []],
    [[-1], [], [12], [], [], [48]],
    [[]],
    [[-12], [], [], []]
]

def test_poly_LC():
    assert poly_LC([]) == 0
    assert poly_LC([2,3,4,5]) == 2
    assert poly_LC([[]]) == []
    assert poly_LC([[2,3,4],[5]]) == [2,3,4]
    assert poly_LC([[[]]]) == [[]]
    assert poly_LC([[[2],[3,4]],[[5]]]) == [[2],[3,4]]

def test_poly_TC():
    assert poly_TC([]) == 0
    assert poly_TC([2,3,4,5]) == 5
    assert poly_TC([[]]) == []
    assert poly_TC([[2,3,4],[5]]) == [5]
    assert poly_TC([[[]]]) == [[]]
    assert poly_TC([[[2],[3,4]],[[5]]]) == [[5]]

def test_poly_nth():
    assert poly_nth([1,2,3], 0) == 3
    assert poly_nth([1,2,3], 1) == 2
    assert poly_nth([1,2,3], 2) == 1

    raises(IndexError, 'poly_nth([3,4,5],-1)')
    raises(IndexError, 'poly_nth([3,4,5], 3)')

def test_poly_level():
    assert poly_level([]) == 1
    assert poly_level([[]]) == 2
    assert poly_level([[[]]]) == 3

    assert poly_level([2,3,4,5]) == 1
    assert poly_level([[2,3,4],[5]]) == 2
    assert poly_level([[[2],[3,4]],[[5]]]) == 3

def test_poly_univariate_p():
    assert poly_univariate_p([]) == True
    assert poly_univariate_p([[]]) == False
    assert poly_univariate_p([[[]]]) == False

    assert poly_univariate_p([2,3,4,5]) == True
    assert poly_univariate_p([[2,3,4],[5]]) == False
    assert poly_univariate_p([[[2],[3,4]],[[5]]]) == False

def test_zzX_values():
    assert zzX_zero(0) == 0
    assert zzX_zero(1) == []
    assert zzX_zero(3) == [[[]]]
    assert zzX_const(1, 7) == [7]
    assert zzX_const(3, 7) == [[[7]]]

    assert zzX_const(0, 3) == 3
    assert zzX_const(3, 0) == [[[]]]

    assert zzX_zero_of(f_0, 1) == [[]]
    assert zzX_zero_of(f_6, 1) == [[[]]]
    assert zzX_const_of(f_0, 7, 1) == [[7]]
    assert zzX_const_of(f_6, 7, 1) == [[[7]]]

    assert zzX_value(3, 1) == [[[1]]]

    assert zzX_value(0, [[1]]) == [[1]]
    assert zzX_value(1, [[1]]) == [[[1]]]
    assert zzX_value(2, [[1]]) == [[[[1]]]]

    assert zzX_zeros_of([1,2,3], 4) == [0,0,0,0]

    assert zzX_zeros_of(f_0, 0, 0) == []
    assert zzX_zeros_of(f_0, 2, 0) == [[[[]]], [[[]]]]
    assert zzX_zeros_of(f_0, 3, 1) == [[[]],[[]],[[]]]

    assert zzX_lift(2, []) == [[[]]]
    assert zzX_lift(2, [[1,2,3], [], [2,3]]) == \
        [[[[1]],[[2]],[[3]]], [[[]]], [[[2]],[[3]]]]

    assert zzX_zero_p([]) == True
    assert zzX_zero_p([[]]) == True
    assert zzX_zero_p([[[]]]) == True
    assert zzX_zero_p([[[1]]]) == False

    assert zzX_one_p([1]) == True
    assert zzX_one_p([[1]]) == True
    assert zzX_one_p([[[1]]]) == True
    assert zzX_one_p([[[12]]]) == False

def test_zzx_degree():
    assert zzx_degree([]) == -1
    assert zzx_degree([1]) == 0
    assert zzx_degree([1,0]) == 1
    assert zzx_degree([1,0,0,0,1]) == 4

def test_zzX_degree():
    assert zzX_degree([[]]) == -1
    assert zzX_degree([[[]]]) == -1

    assert zzX_degree([[1]]) == 0
    assert zzX_degree([[2],[1]]) == 1

    assert zzX_degree(f_0) == 2
    assert zzX_degree(f_1) == 3
    assert zzX_degree(f_2) == 5

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

def test_zzX_strip():
    assert zzX_strip([[]]) == [[]]
    assert zzX_strip([[], []]) == [[]]
    assert zzX_strip([[], [], []]) == [[]]

    assert zzX_strip([[[]]]) == [[[]]]
    assert zzX_strip([[[]], [[]]]) == [[[]]]
    assert zzX_strip([[[]], [[]], [[]]]) == [[[]]]

    assert zzX_strip([[[1]]]) == [[[1]]]
    assert zzX_strip([[[]], [[1]]]) == [[[1]]]
    assert zzX_strip([[[]], [[1]], [[]]]) == [[[1]], [[]]]

def test_zzx_from_to():
    f = [3,0,0,2,0,0,0,0,8]

    g = {8: 3, 5: 2, 0: 8}

    assert zzx_from_dict(g) == f
    assert zzx_to_dict(f) == g

    g = (3*x**8 + 2*x**5 + 8).as_poly(x)

    assert zzx_from_poly(g) == f
    assert zzx_to_poly(f, x) == g

def test_zzX_from_to():
    assert zzX_from_poly(F_0) == f_0
    assert zzX_from_poly(F_1) == f_1
    assert zzX_from_poly(F_2) == f_2
    assert zzX_from_poly(F_3) == f_3
    assert zzX_from_poly(F_4) == f_4
    assert zzX_from_poly(F_5) == f_5
    assert zzX_from_poly(F_6) == f_6

    assert zzX_to_poly(f_0, x,y,z) == F_0
    assert zzX_to_poly(f_1, x,y,z) == F_1
    assert zzX_to_poly(f_2, x,y,z) == F_2
    assert zzX_to_poly(f_3, x,y,z) == F_3
    assert zzX_to_poly(f_4, x,y,z) == F_4
    assert zzX_to_poly(f_5, x,y,z) == F_5
    assert zzX_to_poly(f_6, x,y,z,t) == F_6

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

def test_zzX_add_term():
    assert zzX_add_term(f_0, [[[]]], 3) == f_0

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

def test_zzX_sub_term():
    assert zzX_sub_term(f_0, [[[]]], 3) == f_0

def test_zzx_mul_term():
    assert zzx_mul_term([],    2, 3) == []
    assert zzx_mul_term([1,1], 0, 3) == []

    assert zzx_mul_term([1,2,3], 2, 0) == [2,4,6]
    assert zzx_mul_term([1,2,3], 2, 1) == [2,4,6,0]
    assert zzx_mul_term([1,2,3], 2, 2) == [2,4,6,0,0]
    assert zzx_mul_term([1,2,3], 2, 3) == [2,4,6,0,0,0]

def test_zzX_mul_term():
    assert zzX_mul_term([1,2,3], 2, 1) == [2,4,6,0]

    assert zzX_mul_term([[]], [2], 3) == [[]]
    assert zzX_mul_term([[1]], [], 3) == [[]]

    assert zzX_mul_term([[1,2], [3]], [2], 2) == [[2,4], [6], [], []]

def test_zzx_mul_const():
    assert zzx_mul_const([], 2) == []
    assert zzx_mul_const([1,2,3], 0) == []
    assert zzx_mul_const([1,2,3], 2) == [2,4,6]

def test_zzX_mul_const():
    assert zzX_mul_const(f_0, 2) == [
        [[2,4,6], [4]],
        [[6]],
        [[8,10,12], [2,4,2], [2]]
    ]

def test_zzx_quo_const():
    raises(ZeroDivisionError, 'zzx_quo_const([1,2,3], 0)')
    raises(ExactQuotientFailed, 'zzx_quo_const([1,2,3], 3)')

    assert zzx_quo_const([], 3) == []
    assert zzx_quo_const([6,2,8], 2) == [3,1,4]

def test_zzX_quo_const():
    raises(ExactQuotientFailed, 'zzX_quo_const(f_0, 17)')

def test_zzx_compose_term():
    assert zzx_compose_term([], 17) == []

    assert zzx_compose_term([1,2,3], 1) == [1, 2, 3]
    assert zzx_compose_term([1,2,3], 2) == [1, 0, 2, 0, 3]
    assert zzx_compose_term([1,2,3], 3) == [1, 0, 0, 2, 0, 0, 3]
    assert zzx_compose_term([1,2,3], 4) == [1, 0, 0, 0, 2, 0, 0, 0, 3]

    raises(ValueError, 'zzx_compose_term([1,2,3], 0)')

def test_zzx_abs():
    assert zzx_abs([]) == []
    assert zzx_abs([1]) == [1]
    assert zzx_abs([-7]) == [7]
    assert zzx_abs([-1,2,3]) == [1,2,3]

def test_zzX_abs():
    assert zzX_abs([-1]) == [1]
    assert zzX_abs([[[]]]) == [[[]]]
    assert zzX_abs([[[1]]]) == [[[1]]]
    assert zzX_abs([[[-7]]]) == [[[7]]]

def test_zzx_neg():
    assert zzx_neg([]) == []
    assert zzx_neg([1]) == [-1]
    assert zzx_neg([-7]) == [7]
    assert zzx_neg([-1,2,3]) == [1,-2,-3]

def test_zzX_neg():
    assert zzX_neg([-1]) == [1]
    assert zzX_neg([[[]]]) == [[[]]]
    assert zzX_neg([[[1]]]) == [[[-1]]]
    assert zzX_neg([[[-7]]]) == [[[7]]]

def test_zzx_add():
    assert zzx_add([], []) == []
    assert zzx_add([1], []) == [1]
    assert zzx_add([], [1]) == [1]
    assert zzx_add([1], [1]) == [2]
    assert zzx_add([1], [2]) == [3]

    assert zzx_add([1,2], [1]) == [1,3]
    assert zzx_add([1], [1,2]) == [1,3]

    assert zzx_add([1,2,3], [8,9,10]) == [9,11,13]

def test_zzX_add():
    assert zzX_add([[[]]], [[[]]]) == [[[]]]
    assert zzX_add([[[1]]], [[[]]]) == [[[1]]]
    assert zzX_add([[[]]], [[[1]]]) == [[[1]]]
    assert zzX_add([[[2]]], [[[1]]]) == [[[3]]]
    assert zzX_add([[[1]]], [[[2]]]) == [[[3]]]

def test_zzx_sub():
    assert zzx_sub([], []) == []
    assert zzx_sub([1], []) == [1]
    assert zzx_sub([], [1]) == [-1]
    assert zzx_sub([1], [1]) == []
    assert zzx_sub([1], [2]) == [-1]

    assert zzx_sub([1,2], [1]) == [1,1]
    assert zzx_sub([1], [1,2]) == [-1,-1]

    assert zzx_sub([3,2,1], [8,9,10]) == [-5,-7,-9]

def test_zzX_sub():
    assert zzX_sub([[[]]], [[[]]]) == [[[]]]
    assert zzX_sub([[[1]]], [[[]]]) == [[[1]]]
    assert zzX_sub([[[]]], [[[1]]]) == [[[-1]]]
    assert zzX_sub([[[2]]], [[[1]]]) == [[[1]]]
    assert zzX_sub([[[1]]], [[[2]]]) == [[[-1]]]

def test_zzx_add_sub_mul():
    assert zzx_add_mul([1,2,3], [3,2,1], [1,2]) == [ 3, 9, 7, 5]
    assert zzx_sub_mul([1,2,3], [3,2,1], [1,2]) == [-3,-7,-3, 1]

def test_zzX_add_sub_mul():
    assert zzX_add_mul([[1,2],[3]], [[3],[2,1]], [[1],[2]]) == [[ 3], [ 3,  9], [ 4, 5]]
    assert zzX_sub_mul([[1,2],[3]], [[3],[2,1]], [[1],[2]]) == [[-3], [-1, -5], [-4, 1]]

def test_zzx_mul():
    assert zzx_mul([], []) == []
    assert zzx_mul([], [1]) == []
    assert zzx_mul([1], []) == []
    assert zzx_mul([1], [1]) == [1]
    assert zzx_mul([5], [7]) == [35]

    assert zzx_mul([3,0,0,6,1,2], [4,0,1,0]) == [12,0,3,24,4,14,1,2,0]
    assert zzx_mul([4,0,1,0], [3,0,0,6,1,2]) == [12,0,3,24,4,14,1,2,0]

    assert zzx_mul([2,0,0,1,7], [2,0,0,1,7]) == [4,0,0,4,28,0,1,14,49]

def test_zzX_mul():
    assert zzX_mul([[[]]], [[[]]]) == [[[]]]
    assert zzX_mul([[[1]]], [[[]]]) == [[[]]]
    assert zzX_mul([[[]]], [[[1]]]) == [[[]]]
    assert zzX_mul([[[2]]], [[[1]]]) == [[[2]]]
    assert zzX_mul([[[1]]], [[[2]]]) == [[[2]]]

def test_zzx_sqr():
    assert zzx_sqr([]) == []
    assert zzx_sqr([2]) == [4]
    assert zzx_sqr([1,2]) == [1,4,4]

    assert zzx_sqr([2,0,0,1,7]) == [4,0,0,4,28,0,1,14,49]

def test_zzX_sqr():
    assert zzX_sqr([[[]]]) == [[[]]]
    assert zzX_sqr([[[2]]]) == [[[4]]]

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

def test_zzX_div():
    raises(ZeroDivisionError, "zzX_div([[1,2],[3]], [[]])")
    raises(ZeroDivisionError, "zzX_quo([[1,2],[3]], [[]])")
    raises(ZeroDivisionError, "zzX_rem([[1,2],[3]], [[]])")

    f, g, q, r = [5,4,3,2,1], [1,2,3], [5,-6,0], [20,1]

    assert zzX_div(f, g) == (q, r)
    assert zzX_quo(f, g) == q
    assert zzX_rem(f, g) == r

    f, g, q, r = [[[1]]], f_0, [[[]]], [[[1]]]

    assert zzX_div(f, g) == (q, r)
    assert zzX_quo(f, g) == q
    assert zzX_rem(f, g) == r

def test_zzx_gcd():
    assert zzx_heu_gcd([], []) == ([], [], [])
    assert zzx_mod_gcd([], []) == ([], [], [])

    assert zzx_heu_gcd([2], []) == ([2], [1], [])
    assert zzx_mod_gcd([2], []) == ([2], [1], [])

    assert zzx_heu_gcd([], [2,4]) == ([2,4], [], [1])
    assert zzx_mod_gcd([], [2,4]) == ([2,4], [], [1])

    assert zzx_heu_gcd([2,4], []) == ([2,4], [1], [])
    assert zzx_mod_gcd([2,4], []) == ([2,4], [1], [])

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

    assert zzx_gcd(f, g) == h
    assert zzx_cofactors(f, g) == (h, cff, cfg)

def test_zzx_heu_gcd():
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

    assert zzx_heu_gcd(f, zzx_diff(f))[0] == g

def test_zzX_heu_gcd():
    f,g,h = zzX_fateman_poly_F_1(2)
    H, cff, cfg = zzX_heu_gcd(f, g)

    assert H == h and zzX_mul(H, cff) == f \
                  and zzX_mul(H, cfg) == g

    f,g,h = zzX_fateman_poly_F_1(4)
    H, cff, cfg = zzX_heu_gcd(f, g)

    assert H == h and zzX_mul(H, cff) == f \
                  and zzX_mul(H, cfg) == g

    f,g,h = zzX_fateman_poly_F_1(6)
    H, cff, cfg = zzX_heu_gcd(f, g)

    assert H == h and zzX_mul(H, cff) == f \
                  and zzX_mul(H, cfg) == g

    f,g,h = zzX_fateman_poly_F_2(2)
    H, cff, cfg = zzX_heu_gcd(f, g)

    assert H == h and zzX_mul(H, cff) == f \
                  and zzX_mul(H, cfg) == g

    f,g,h = zzX_fateman_poly_F_3(2)
    H, cff, cfg = zzX_heu_gcd(f, g)

    assert H == h and zzX_mul(H, cff) == f \
                  and zzX_mul(H, cfg) == g

    f,g,h = zzX_fateman_poly_F_1(8)

    try:
        import gmpy

        H, cff, cfg = zzX_heu_gcd(f, g)

        assert H == h and zzX_mul(H, cff) == f \
                      and zzX_mul(H, cfg) == g
    except ImportError:
        raises(HeuristicGCDFailed, 'zzX_heu_gcd(f, g)')

def test_zzx_norm():
    assert zzx_max_norm([]) == 0
    assert zzx_max_norm([1]) == 1
    assert zzx_max_norm([1,4,2,3]) == 4

    assert zzx_l1_norm([]) == 0
    assert zzx_l1_norm([1]) == 1
    assert zzx_l1_norm([1,4,2,3]) == 10

def test_zzX_norm():
    assert zzX_max_norm([[[]]]) == 0
    assert zzX_max_norm([[[1]]]) == 1

    assert zzX_l1_norm([[[]]]) == 0
    assert zzX_l1_norm([[[1]]]) == 1

    assert zzX_max_norm(f_0) == 6
    assert zzX_l1_norm(f_0) == 31

def test_zzx_diff():
    assert zzx_diff([]) == []
    assert zzx_diff([7]) == []
    assert zzx_diff([2,7]) == [2]
    assert zzx_diff([1,2,1]) == [2,2]
    assert zzx_diff([1,2,3,4]) == [3,4,3]
    assert zzx_diff([1,-1,0,0,2]) == [4,-3,0,0]

def test_zzX_diff():
    assert zzX_diff([]) == []

    assert zzX_diff([[]]) == [[]]

    assert zzX_diff([[[]]]) == [[[]]]
    assert zzX_diff([[[1], [2]]]) == [[[]]]

    assert zzX_diff([[[1]], [[]]]) == [[[1]]]
    assert zzX_diff([[[3]], [[1]], [[]]]) == [[[6]], [[1]]]

def test_zzx_eval():
    assert zzx_eval([], 7) == 0
    assert zzx_eval([1], 0) == 1
    assert zzx_eval([1,2,3], 7) == 66

def test_zzX_eval():
    assert zzX_eval([], 3) == 0

    assert zzX_eval([[]], 3) == []
    assert zzX_eval([[[]]], 3) == [[]]

    assert zzX_eval([[1,2]], 0) == [1,2]

    assert zzX_eval([[[1]]], 3) == [[1]]
    assert zzX_eval([[[1, 2]]], 3) == [[1, 2]]

    assert zzX_eval([[3, 2], [1, 2]], 3) == [10, 8]
    assert zzX_eval([[[3, 2]], [[1, 2]]], 3) == [[10, 8]]

    assert zzX_eval_coeffs([], []) == []
    assert zzX_eval_coeffs([[]], []) == []
    assert zzX_eval_coeffs([[[]]], []) == []

    assert zzX_eval_coeffs(f_0, [-17, 8]) == [-1409, 3, 85902]
    assert zzX_eval_coeffs(f_1, [-17, 8]) == [-136, 15699, 9166, -27144]

    assert zzX_eval_coeffs(f_2, [-12, 3]) == [-1377, 0, -702, -1224, 0, -624]
    assert zzX_eval_coeffs(f_3, [-12, 3]) == [144, 82, -5181, -28872, -14868, -540]

    assert zzX_eval_coeffs(f_4, [25, -1]) == [152587890625, 9765625, -59605407714843750,
        -3839159765625, -1562475, 9536712644531250, 610349546750, -4, 24414375000, 1562520]
    assert zzX_eval_coeffs(f_5, [25, -1]) == [-1, -78, -2028, -17576]

    assert zzX_eval_coeffs(f_6, [0, 2, 4]) == [5040, 0, 0, 4480]

def test_zzx_trunc():
    assert zzx_trunc([1,2,3,4,5,6], 3) == [1, -1, 0, 1, -1, 0]
    assert zzx_trunc([6,5,4,3,2,1], 3) ==    [-1, 1, 0, -1, 1]

def test_zzX_trunc():
    assert zzX_trunc(f_0, 3) == \
        [[[1, -1, 0], [-1]], [[]], [[1, -1, 0], [1, -1, 1], [1]]]

def test_zzX_zz_LC_TC():
    assert zzX_zz_LC(f_0) == 1
    assert zzX_zz_TC(f_0) == 1

    assert zzX_zz_LC(f_1) == 1
    assert zzX_zz_TC(f_1) == 6000

    assert zzX_zz_LC(f_2) == 1
    assert zzX_zz_TC(f_2) == -990

    assert zzX_zz_LC(f_3) == 1
    assert zzX_zz_TC(f_3) == 0

    assert zzX_zz_LC(f_4) == -1
    assert zzX_zz_TC(f_4) == 0

    assert zzX_zz_LC(f_5) == -1
    assert zzX_zz_TC(f_5) == 0

    assert zzX_zz_LC(f_6) == 2115
    assert zzX_zz_TC(f_6) == 0

def test_zzx_content():
    assert zzx_content([]) == 0
    assert zzx_content([1]) == 1
    assert zzx_content([1,1]) == 1
    assert zzx_content([2,2]) == 2
    assert zzx_content([1,2,1]) == 1
    assert zzx_content([2,4,2]) == 2

def test_zzX_content():
    f, g, F = [3,2,1], [1], []

    for i in xrange(0, 5):
        g = zzX_mul(g, f)
        F.insert(0, g)

    assert zzX_content(F) == f

    assert zzX_one_p(zzX_content(f_4))
    assert zzX_one_p(zzX_content(f_5))
    assert zzX_one_p(zzX_content(f_6))

def test_zzX_zz_content():
    assert zzX_zz_content(f_0) == 1
    assert zzX_zz_content(zzX_mul_const(f_0, 2)) == 2

    assert zzX_zz_content(f_1) == 1
    assert zzX_zz_content(zzX_mul_const(f_1, 3)) == 3

    assert zzX_zz_content(f_2) == 1
    assert zzX_zz_content(zzX_mul_const(f_2, 4)) == 4

    assert zzX_zz_content(f_3) == 1
    assert zzX_zz_content(zzX_mul_const(f_3, 5)) == 5

    assert zzX_zz_content(f_4) == 1
    assert zzX_zz_content(zzX_mul_const(f_4, 6)) == 6

    assert zzX_zz_content(f_5) == 1
    assert zzX_zz_content(zzX_mul_const(f_5, 7)) == 7

    assert zzX_zz_content(f_6) == 1
    assert zzX_zz_content(zzX_mul_const(f_6, 8)) == 8

def test_zzx_primitive():
    assert zzx_primitive([]) == (0, [])
    assert zzx_primitive([1]) == (1, [1])
    assert zzx_primitive([1,1]) == (1, [1,1])
    assert zzx_primitive([2,2]) == (2, [1,1])
    assert zzx_primitive([1,2,1]) == (1, [1,2,1])
    assert zzx_primitive([2,4,2]) == (2, [1,2,1])

def test_zzX_primitive():
    f, g, F = [3,2,1], [1], []

    for i in xrange(0, 5):
        g = zzX_mul(g, f)
        F.insert(0, g)

    assert zzX_primitive(F) == (f,
        [ zzX_quo(cf, f) for cf in F ])

    cont, f = zzX_primitive(f_4)
    assert zzX_one_p(cont) and f == f_4
    cont, f = zzX_primitive(f_5)
    assert zzX_one_p(cont) and f == f_5
    cont, f = zzX_primitive(f_6)
    assert zzX_one_p(cont) and f == f_6

def test_zzX_zz_primitive():
    assert zzX_zz_primitive(f_0) == (1, f_0)
    assert zzX_zz_primitive(zzX_mul_const(f_0, 2)) == (2, f_0)

    assert zzX_zz_primitive(f_1) == (1, f_1)
    assert zzX_zz_primitive(zzX_mul_const(f_1, 3)) == (3, f_1)

    assert zzX_zz_primitive(f_2) == (1, f_2)
    assert zzX_zz_primitive(zzX_mul_const(f_2, 4)) == (4, f_2)

    assert zzX_zz_primitive(f_3) == (1, f_3)
    assert zzX_zz_primitive(zzX_mul_const(f_3, 5)) == (5, f_3)

    assert zzX_zz_primitive(f_4) == (1, f_4)
    assert zzX_zz_primitive(zzX_mul_const(f_4, 6)) == (6, f_4)

    assert zzX_zz_primitive(f_5) == (1, f_5)
    assert zzX_zz_primitive(zzX_mul_const(f_5, 7)) == (7, f_5)

    assert zzX_zz_primitive(f_6) == (1, f_6)
    assert zzX_zz_primitive(zzX_mul_const(f_6, 8)) == (8, f_6)

def test_zzx_sqf():
    assert zzx_sqf_part([2,2]) == [1,1]
    assert zzx_sqf_p([2,2]) == True

    assert zzx_sqf_part([1,0,1,1]) == [1,0,1,1]
    assert zzx_sqf_p([1,0,1,1]) == True

    assert zzx_sqf_part([2,3,0,0]) == [2,3,0]
    assert zzx_sqf_p([2,3,0,0]) == False

def test_zzX_sqf():
    assert zzX_sqf_p(f_0) == True
    assert zzX_sqf_p(zzX_sqr(f_0)) == False
    assert zzX_sqf_p(f_1) == True
    assert zzX_sqf_p(zzX_sqr(f_1)) == False
    assert zzX_sqf_p(f_2) == True
    assert zzX_sqf_p(zzX_sqr(f_2)) == False
    assert zzX_sqf_p(f_3) == True
    assert zzX_sqf_p(zzX_sqr(f_3)) == False
    assert zzX_sqf_p(f_5) == False
    assert zzX_sqf_p(zzX_sqr(f_5)) == False

    assert zzX_sqf_p(f_4) == True
    assert zzX_sqf_part(f_4) == f_4
    assert zzX_sqf_p(f_6) == True
    assert zzX_sqf_part(f_6) == f_6

    assert zzX_sqf_part(f_5) == [[[-1]], [[-1], [1, 0]]]

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
    assert zzx_factor([ ]) == (0, [])
    assert zzx_factor([7]) == (7, [])

    f = [1,0,0,1,1]

    for i in xrange(0, 20):
        assert zzx_factor(f) == (1, [(f, 1)])

    assert zzx_factor([2,4]) == \
        (2, [([1, 2], 1)])

    assert zzx_factor([18,12,2]) == \
        (2, [([3, 1], 2)])

    assert zzx_factor([-9,0,1]) == \
        (-1, [([3,-1], 1),
              ([3, 1], 1)])

    assert zzx_factor([1,-6,11,-6]) == \
        (1, [([1,-3], 1),
             ([1,-2], 1),
             ([1,-1], 1)])

    assert zzx_factor([-1,0,0,0,1,0,0]) == \
        (-1, [([1,-1], 1),
              ([1, 1], 1),
              ([1, 0], 2),
              ([1, 0, 1], 1)])

    f = [1080, 5184, 2099, 744, 2736, -648, 129, 0, -324]

    assert zzx_factor(f) == \
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

    assert zzx_factor(f) == \
        (-95367431640625, [([5, -1], 1),
                           ([100, 10, -1], 2),
                           ([625, 125, 25, 5, 1], 1),
                           ([10000, -3000, 400, -20, 1], 2),
                           ([10000,  2000, 400,  30, 1], 2)])

    f = zzx_from_dict({10:1, 0:-1})

    F_0 = zzx_factor(f, cyclotomic=True)
    F_1 = zzx_factor(f, cyclotomic=False)

    assert F_0 == F_1 == \
        (1, [([1,-1], 1),
             ([1, 1], 1),
             ([1,-1, 1,-1, 1], 1),
             ([1, 1, 1, 1, 1], 1)])

    f = zzx_from_dict({10:1, 0:1})

    F_0 = zzx_factor(f, cyclotomic=True)
    F_1 = zzx_factor(f, cyclotomic=False)

    assert F_0 == F_1 == \
        (1, [([1, 0, 1], 1),
             ([1, 0, -1, 0, 1, 0, -1, 0, 1], 1)])

def test_zzx_cyclotomic_factor():
    assert zzx_cyclotomic_factor([]) is None
    assert zzx_cyclotomic_factor([1]) is None

    f = zzx_from_dict({10:2, 0:-1})
    assert zzx_cyclotomic_factor(f) is None
    f = zzx_from_dict({10:1, 0:-3})
    assert zzx_cyclotomic_factor(f) is None
    f = zzx_from_dict({10:1, 5:1, 0:-1})
    assert zzx_cyclotomic_factor(f) is None

    f = zzx_from_dict({1:1,0:1})
    assert zzx_cyclotomic_factor(f) == \
         [[1, 1]]

    f = zzx_from_dict({1:1,0:-1})
    assert zzx_cyclotomic_factor(f) == \
        [[1, -1]]

    f = zzx_from_dict({2:1,0:1})
    assert zzx_cyclotomic_factor(f) == \
        [[1, 0, 1]]

    f = zzx_from_dict({2:1,0:-1})
    assert zzx_cyclotomic_factor(f) == \
        [[1,-1],
         [1, 1]]

    f = zzx_from_dict({27:1,0:1})
    assert zzx_cyclotomic_factor(f) == \
        [[1, 1],
         [1, -1, 1],
         [1, 0, 0, -1, 0, 0, 1],
         [1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1]]

    f = zzx_from_dict({27:1,0:-1})
    assert zzx_cyclotomic_factor(f) == \
        [[1, -1],
         [1, 1, 1],
         [1, 0, 0, 1, 0, 0, 1],
         [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]]

