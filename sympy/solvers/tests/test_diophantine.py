from sympy.solvers.diophantine import (diop_solve, diop_DN, diop_bf_DN, length, transformation_to_DN, find_DN, equivalent,
    square_factor, pairwise_prime, descent,
    ldescent, diophantine, transformation_to_normal, sum_of_four_squares, sum_of_three_squares,
    prime_as_sum_of_two_squares, partition, power_representation)

from sympy import symbols, Integer, Matrix, simplify, Subs, S, factor_list
from sympy.core.function import _mexpand
from sympy.core.compatibility import range
from sympy.functions.elementary.trigonometric import sin
from sympy.utilities.pytest import slow, raises
from sympy.utilities import default_sort_key

x, y, z, w, t, X, Y, Z = symbols("x, y, z, w, t, X, Y, Z", integer=True)

def test_input_format():
    raises(TypeError, lambda: diophantine(sin(x)))

def test_univariate():

    assert diop_solve((x - 1)*(x - 2)**2) == set([(Integer(1),), (Integer(2),)])
    assert diop_solve((x - 1)*(x - 2)) == set([(Integer(1),), (Integer(2),)])


def test_linear():

    assert diop_solve(2*x + 3*y - 5) == (3*t - 5, -2*t + 5)
    assert diop_solve(3*y + 2*x - 5) == (3*t - 5, -2*t + 5)
    assert diop_solve(2*x - 3*y - 5) == (-3*t - 5, -2*t - 5)
    assert diop_solve(-2*x - 3*y - 5) == (-3*t + 5, 2*t - 5)
    assert diop_solve(7*x + 5*y) == (5*t, -7*t)
    assert diop_solve(2*x + 4*y) == (2*t, -t)
    assert diop_solve(4*x + 6*y - 4) == (3*t - 2, -2*t + 2)
    assert diop_solve(4*x + 6*y - 3) == (None, None)
    assert diop_solve(4*x + 3*y -4*z + 5) == \
           (3*t + 4*z - 5, -4*t - 4*z + 5, z)
    assert diop_solve(4*x + 2*y + 8*z - 5) == (None, None, None)
    assert diop_solve(5*x + 7*y - 2*z - 6) == \
           (7*t + 6*z + 18, -5*t - 4*z - 12, z)
    assert diop_solve(3*x - 6*y + 12*z - 9) == \
           (-2*t - 4*z + 3, -t, z)
    assert diop_solve(x + 3*y - 4*z + w - 6) == \
           (t, -t - 3*y + 4*z + 6, y, z)


def test_quadratic_simple_hyperbolic_case():

    # Simple Hyperbolic case: A = C = 0 and B != 0
    assert diop_solve(3*x*y + 34*x - 12*y + 1) == \
        set([(-Integer(133), -Integer(11)), (Integer(5), -Integer(57))])
    assert diop_solve(6*x*y + 2*x + 3*y + 1) == set([])
    assert diop_solve(-13*x*y + 2*x - 4*y - 54) == set([(Integer(27), Integer(0))])
    assert diop_solve(-27*x*y - 30*x - 12*y - 54) == set([(-Integer(14), -Integer(1))])
    assert diop_solve(2*x*y + 5*x + 56*y + 7) == set([(-Integer(161), -Integer(3)),\
        (-Integer(47),-Integer(6)), (-Integer(35), -Integer(12)), (-Integer(29), -Integer(69)),\
        (-Integer(27), Integer(64)), (-Integer(21), Integer(7)),(-Integer(9), Integer(1)),\
        (Integer(105), -Integer(2))])
    assert diop_solve(6*x*y + 9*x + 2*y + 3) == set([])
    assert diop_solve(x*y + x + y + 1) == set([(-Integer(1), t), (t, -Integer(1))])
    assert diophantine(48*x*y)


def test_quadratic_elliptical_case():

    # Elliptical case: B**2 - 4AC < 0
    # Two test cases highlighted require lot of memory due to quadratic_congruence() method.
    # This above method should be replaced by Pernici's square_mod() method when his PR gets merged.

    #assert diop_solve(42*x**2 + 8*x*y + 15*y**2 + 23*x + 17*y - 4915) == set([(-Integer(11), -Integer(1))])
    assert diop_solve(4*x**2 + 3*y**2 + 5*x - 11*y + 12) == set([])
    assert diop_solve(x**2 + y**2 + 2*x + 2*y + 2) == set([(-Integer(1), -Integer(1))])
    #assert diop_solve(15*x**2 - 9*x*y + 14*y**2 - 23*x - 14*y - 4950) == set([(-Integer(15), Integer(6))])
    assert diop_solve(10*x**2 + 12*x*y + 12*y**2 - 34) == \
        set([(Integer(1), -Integer(2)), (-Integer(1), -Integer(1)),(Integer(1), Integer(1)), (-Integer(1), Integer(2))])


def test_quadratic_parabolic_case():

    # Parabolic case: B**2 - 4AC = 0
    assert diop_solve(8*x**2 - 24*x*y + 18*y**2 + 5*x + 7*y + 16) == \
        set([(-174*t**2 + 17*t - 2, -116*t**2 + 21*t - 2), (-174*t**2 + 41*t - 4, -116*t**2 + 37*t - 4)])
    assert diop_solve(8*x**2 - 24*x*y + 18*y**2 + 6*x + 12*y - 6) == \
        set([(-63*t**2 + 12*t, -42*t**2 + 15*t -1), (-63*t**2 + 30*t - 3, -42*t**2 + 27*t - 4)])
    assert diop_solve(8*x**2 + 24*x*y + 18*y**2 + 4*x + 6*y - 7) == set([])
    assert diop_solve(x**2 + 2*x*y + y**2 + 2*x + 2*y + 1) == set([(t,-t - 1)])
    assert diop_solve(x**2 - 2*x*y + y**2 + 2*x + 2*y + 1) == \
        set([(-4*t**2, -4*t**2 + 4*t - 1),(-4*t**2 + 4*t -1, -4*t**2 + 8*t - 4)])
    assert check_solutions(y**2 - 41*x + 40)


def test_quadratic_perfect_square():

    # B**2 - 4*A*C > 0
    # B**2 - 4*A*C is a perfect square
    assert diop_solve(48*x*y) == set([(Integer(0), t), (t, Integer(0))])
    assert diop_solve(4*x**2 - 5*x*y + y**2 + 2) == \
        set([(-Integer(1), -Integer(3)),(-Integer(1), -Integer(2)),(Integer(1), Integer(2)),(Integer(1), Integer(3))])
    assert diop_solve(-2*x**2 - 3*x*y + 2*y**2 -2*x - 17*y + 25) == set([(Integer(4), Integer(15))])
    assert diop_solve(12*x**2 + 13*x*y + 3*y**2 - 2*x + 3*y - 12) == \
        set([(-Integer(6), Integer(9)), (-Integer(2), Integer(5)), (Integer(4), -Integer(4)), (-Integer(6), Integer(16))])
    assert diop_solve(8*x**2 + 10*x*y + 2*y**2 - 32*x - 13*y - 23) == \
        set([(-Integer(44), Integer(47)), (Integer(22), -Integer(85))])
    assert diop_solve(4*x**2 - 4*x*y - 3*y- 8*x - 3) == \
        set([(-Integer(1), -Integer(9)), (-Integer(6), -Integer(9)), (Integer(0), -Integer(1)), (Integer(1), -Integer(1))])
    assert diop_solve(- 4*x*y - 4*y**2 - 3*y- 5*x - 10) == \
        set([(-Integer(2), Integer(0)), (-Integer(11), -Integer(1)), (-Integer(5), Integer(5))])
    assert diop_solve(x**2 - y**2 - 2*x - 2*y) == set([(t, -t), (-t, -t - 2)])
    assert diop_solve(x**2 - 9*y**2 - 2*x - 6*y) == set([(-3*t + 2, -t), (3*t, -t)])
    assert diop_solve(4*x**2 - 9*y**2 - 4*x - 12*y - 3) == set([(-3*t - 3, -2*t - 3), (3*t + 1, -2*t - 1)])


def test_quadratic_non_perfect_square():

    # B**2 - 4*A*C is not a perfect square
    # Used check_solutions() since the solutions are complex expressions involving
    # square roots and exponents
    assert check_solutions(x**2 - 2*x - 5*y**2)
    assert check_solutions(3*x**2 - 2*y**2 - 2*x - 2*y)
    assert check_solutions(x**2 - x*y - y**2 - 3*y)
    assert check_solutions(x**2 - 9*y**2 - 2*x - 6*y)

def test_issue_9106():
    assert check_integrality(-48 - 2*x*(3*x - 1) + y*(3*y - 1))

@slow
def test_quadratic_non_perfect_slow():

    assert check_solutions(8*x**2 + 10*x*y - 2*y**2 - 32*x - 13*y - 23)
    # This leads to very large numbers.
    # assert check_solutions(5*x**2 - 13*x*y + y**2 - 4*x - 4*y - 15)
    assert check_solutions(-3*x**2 - 2*x*y + 7*y**2 - 5*x - 7)
    assert check_solutions(-4 - x + 4*x**2 - y - 3*x*y - 4*y**2)
    assert check_solutions(1 + 2*x + 2*x**2 + 2*y + x*y - 2*y**2)


def test_DN():

    # Most of the test cases were adapted from,
    # Solving the generalized Pell equation x**2 - D*y**2 = N, John P. Robertson, July 31, 2004.
    # http://www.jpr2718.org/pell.pdf
    # others are verified using Wolfram Alpha.

    # Covers cases where D <= 0 or D > 0 and D is a square or N = 0
    # Solutions are straightforward in these cases.
    assert diop_DN(3, 0) == [(0, 0)]
    assert diop_DN(-17, -5) == []
    assert diop_DN(-19, 23) == [(2, 1)]
    assert diop_DN(-13, 17) == [(2, 1)]
    assert diop_DN(-15, 13) == []
    assert diop_DN(0, 5) == []
    assert diop_DN(0, 9) == [(3, t)]
    assert diop_DN(9, 0) == [(3*t, t)]
    assert diop_DN(16, 24) == []
    assert diop_DN(9, 180) == [(18, 4)]
    assert diop_DN(9, -180) == [(12, 6)]
    assert diop_DN(7, 0) == [(0, 0)]

    # When equation is x**2 + y**2 = N
    # Solutions are interchangeable
    assert diop_DN(-1, 5) == [(2, 1)]
    assert diop_DN(-1, 169) == [(12, 5), (0, 13)]

    # D > 0 and D is not a square

    # N = 1
    assert diop_DN(13, 1) == [(649, 180)]
    assert diop_DN(980, 1) == [(51841, 1656)]
    assert diop_DN(981, 1) == [(158070671986249, 5046808151700)]
    assert diop_DN(986, 1) == [(49299, 1570)]
    assert diop_DN(991, 1) == [(379516400906811930638014896080, 12055735790331359447442538767)]
    assert diop_DN(17, 1) == [(33, 8)]
    assert diop_DN(19, 1) == [(170, 39)]

    # N = -1
    assert diop_DN(13, -1) == [(18, 5)]
    assert diop_DN(991, -1) == []
    assert diop_DN(41, -1) == [(32, 5)]
    assert diop_DN(290, -1) == [(17, 1)]
    assert diop_DN(21257, -1) == [(13913102721304, 95427381109)]
    assert diop_DN(32, -1) == []

    # |N| > 1
    # Some tests were created using calculator at
    # http://www.numbertheory.org/php/patz.html

    assert diop_DN(13, -4) == [(3, 1), (393, 109), (36, 10)]
    # Source I referred returned (3, 1), (393, 109) and (-3, 1) as fundamental solutions
    # So (-3, 1) and (393, 109) should be in the same equivalent class
    assert equivalent(-3, 1, 393, 109, 13, -4) == True

    assert diop_DN(13, 27) == [(220, 61), (40, 11), (768, 213), (12, 3)]
    assert set(diop_DN(157, 12)) == \
    set([(Integer(13), Integer(1)), (Integer(10663), Integer(851)), (Integer(579160), Integer(46222)), \
        (Integer(483790960),Integer(38610722)), (Integer(26277068347), Integer(2097138361)), (Integer(21950079635497), Integer(1751807067011))])
    assert diop_DN(13, 25) == [(3245, 900)]
    assert diop_DN(192, 18) == []
    assert diop_DN(23, 13) == [(-6, 1), (6, 1)]
    assert diop_DN(167, 2) == [(13, 1)]
    assert diop_DN(167, -2) == []

    assert diop_DN(123, -2) == [(11, 1)]
    # One calculator returned [(11, 1), (-11, 1)] but both of these are in
    # the same equivalence class
    assert equivalent(11, 1, -11, 1, 123, -2)

    assert diop_DN(123, -23) == [(-10, 1), (10, 1)]


def test_bf_pell():

    assert diop_bf_DN(13, -4) == [(3, 1), (-3, 1), (36, 10)]
    assert diop_bf_DN(13, 27) == [(12, 3), (-12, 3), (40, 11), (-40, 11)]
    assert diop_bf_DN(167, -2) == []
    assert diop_bf_DN(1729, 1) == [(44611924489705, 1072885712316)]
    assert diop_bf_DN(89, -8) == [(9, 1), (-9, 1)]
    assert diop_bf_DN(21257, -1) == [(13913102721304, 95427381109)]
    assert diop_bf_DN(340, -4) == [(756, 41)]


def test_length():

    assert length(-2, 4, 5) == 3
    assert length(-5, 4, 17) == 4
    assert length(0, 4, 13) == 6
    assert length(-31, 8, 613) == 67
    assert length(7, 13, 11) == 23
    assert length(-40, 5, 23) == 4


def is_pell_transformation_ok(eq):
    """
    Test whether X*Y, X, or Y terms are present in the equation
    after transforming the equation using the transformation returned
    by transformation_to_pell(). If they are not present we are good.
    Moreover, coefficient of X**2 should be a divisor of coefficient of
    Y**2 and the constant term.
    """
    A, B = transformation_to_DN(eq)
    u = (A*Matrix([X, Y]) + B)[0]
    v = (A*Matrix([X, Y]) + B)[1]
    simplified = _mexpand(Subs(eq, (x, y), (u, v)).doit())

    coeff = dict([reversed(t.as_independent(*[X, Y])) for t in simplified.args])

    for term in [X*Y, X, Y]:
        if term in coeff.keys():
            return False

    for term in [X**2, Y**2, Integer(1)]:
        if term not in coeff.keys():
            coeff[term] = Integer(0)

    if coeff[X**2] != 0:
        return isinstance(S(coeff[Y**2])/coeff[X**2], Integer) and isinstance(S(coeff[Integer(1)])/coeff[X**2], Integer)

    return True


def test_transformation_to_pell():

    assert is_pell_transformation_ok(-13*x**2 - 7*x*y + y**2 + 2*x - 2*y - 14)
    assert is_pell_transformation_ok(-17*x**2 + 19*x*y - 7*y**2 - 5*x - 13*y - 23)
    assert is_pell_transformation_ok(x**2 - y**2 + 17)
    assert is_pell_transformation_ok(-x**2 + 7*y**2 - 23)
    assert is_pell_transformation_ok(25*x**2 - 45*x*y + 5*y**2 - 5*x - 10*y + 5)
    assert is_pell_transformation_ok(190*x**2 + 30*x*y + y**2 - 3*y - 170*x - 130)
    assert is_pell_transformation_ok(x**2 - 2*x*y -190*y**2 - 7*y - 23*x - 89)
    assert is_pell_transformation_ok(15*x**2 - 9*x*y + 14*y**2 - 23*x - 14*y - 4950)


def test_find_DN():

    assert find_DN(x**2 - 2*x - y**2) == (1, 1)
    assert find_DN(x**2 - 3*y**2 - 5) == (3, 5)
    assert find_DN(x**2 - 2*x*y - 4*y**2 - 7) == (5, 7)
    assert find_DN(4*x**2 - 8*x*y - y**2 - 9) == (20, 36)
    assert find_DN(7*x**2 - 2*x*y - y**2 - 12) == (8, 84)
    assert find_DN(-3*x**2 + 4*x*y -y**2) == (1, 0)
    assert find_DN(-13*x**2 - 7*x*y + y**2 + 2*x - 2*y -14) == (101, -7825480)


def test_ldescent():

    # Equations which have solutions
    u = ([(13, 23), (3, -11), (41, -113), (4, -7), (-7, 4), (91, -3), (1, 1), (1, -1),
        (4, 32), (17, 13), (123689, 1), (19, -570)])
    for a, b in u:
        w, x, y = ldescent(a, b)
        assert a*x**2 + b*y**2 == w**2


def test_diop_ternary_quadratic_normal():

    assert check_solutions(234*x**2 - 65601*y**2 - z**2)
    assert check_solutions(23*x**2 + 616*y**2 - z**2)
    assert check_solutions(5*x**2 + 4*y**2 - z**2)
    assert check_solutions(3*x**2 + 6*y**2 - 3*z**2)
    assert check_solutions(x**2 + 3*y**2 - z**2)
    assert check_solutions(4*x**2 + 5*y**2 - z**2)
    assert check_solutions(x**2 + y**2 - z**2)
    assert check_solutions(16*x**2 + y**2 - 25*z**2)
    assert check_solutions(6*x**2 - y**2 + 10*z**2)
    assert check_solutions(213*x**2 + 12*y**2 - 9*z**2)
    assert check_solutions(34*x**2 - 3*y**2 - 301*z**2)
    assert check_solutions(124*x**2 - 30*y**2 - 7729*z**2)


def is_normal_transformation_ok(eq):

    A = transformation_to_normal(eq)
    X, Y, Z = A*Matrix([x, y, z])
    simplified = _mexpand(Subs(eq, (x, y, z), (X, Y, Z)).doit())

    coeff = dict([reversed(t.as_independent(*[X, Y, Z])) for t in simplified.args])
    for term in [X*Y, Y*Z, X*Z]:
        if term in coeff.keys():
            return False

    return True


def test_transformation_to_normal():

    assert is_normal_transformation_ok(x**2 + 3*y**2 + z**2 - 13*x*y - 16*y*z + 12*x*z)
    assert is_normal_transformation_ok(x**2 + 3*y**2 - 100*z**2)
    assert is_normal_transformation_ok(x**2 + 23*y*z)
    assert is_normal_transformation_ok(3*y**2 - 100*z**2 - 12*x*y)
    assert is_normal_transformation_ok(x**2 + 23*x*y - 34*y*z + 12*x*z)
    assert is_normal_transformation_ok(z**2 + 34*x*y - 23*y*z + x*z)
    assert is_normal_transformation_ok(x**2 + y**2 + z**2 - x*y - y*z - x*z)


def test_diop_ternary_quadratic():
    # Commented out test cases should be uncommented after
    # the bug with factor_list() gets merged.

    assert check_solutions(2*x**2 + z**2 + y**2 - 4*x*y)
    assert check_solutions(x**2 - y**2 - z**2 - x*y - y*z)
    assert check_solutions(3*x**2 - x*y - y*z - x*z)
    assert check_solutions(x**2 - y*z - x*z)
    #assert check_solutions(5*x**2 - 3*x*y - x*z)
    assert check_solutions(4*x**2 - 5*y**2 - x*z)
    assert check_solutions(3*x**2 + 2*y**2 - z**2 - 2*x*y + 5*y*z - 7*y*z)
    assert check_solutions(8*x**2 - 12*y*z)
    assert check_solutions(45*x**2 - 7*y**2 - 8*x*y - z**2)
    assert check_solutions(x**2 - 49*y**2 - z**2 + 13*z*y -8*x*y)
    assert check_solutions(90*x**2 + 3*y**2 + 5*x*y + 2*z*y + 5*x*z)
    assert check_solutions(x**2 + 3*y**2 + z**2 - x*y - 17*y*z)
    assert check_solutions(x**2 + 3*y**2 + z**2 - x*y - 16*y*z + 12*x*z)
    assert check_solutions(x**2 + 3*y**2 + z**2 - 13*x*y - 16*y*z + 12*x*z)
    assert check_solutions(x*y - 7*y*z + 13*x*z)


def test_pairwise_prime():

    assert pairwise_prime(6, 10, 15) == (5, 3, 2)
    assert pairwise_prime(2, 3, 5) == (2, 3, 5)
    assert pairwise_prime(1, 4, 7) == (1, 4, 7)
    assert pairwise_prime(4, 6, 5) == (1, 6, 5)
    assert pairwise_prime(6, 10, -15) == (5, 3, -2)
    assert pairwise_prime(-6, -10, -15) == (-5, -3, -2)
    assert pairwise_prime(4, -6, -5) == (1, -6, -5)


def test_square_factor():

    assert square_factor(1) == square_factor(-1) == 1
    assert square_factor(0) == 1
    assert square_factor(5) == square_factor(-5) == 1
    assert square_factor(4) == square_factor(-4) == 2
    assert square_factor(12) == square_factor(-12) == 2
    assert square_factor(6) == 1
    assert square_factor(18) == 3
    assert square_factor(52) == 2
    assert square_factor(49) == 7
    assert square_factor(392) == 14


def test_parametrize_ternary_quadratic():

    assert check_solutions(x**2 + y**2 - z**2)
    assert check_solutions(x**2 + 2*x*y + z**2)
    assert check_solutions(234*x**2 - 65601*y**2 - z**2)
    assert check_solutions(3*x**2 + 2*y**2 - z**2 - 2*x*y + 5*y*z - 7*y*z)
    assert check_solutions(x**2 - y**2 - z**2)
    assert check_solutions(x**2 - 49*y**2 - z**2 + 13*z*y - 8*x*y)
    assert check_solutions(8*x*y + z**2)
    assert check_solutions(124*x**2 - 30*y**2 - 7729*z**2)
    assert check_solutions(236*x**2 - 225*y**2 - 11*x*y - 13*y*z - 17*x*z)
    assert check_solutions(90*x**2 + 3*y**2 + 5*x*y + 2*z*y + 5*x*z)
    assert check_solutions(124*x**2 - 30*y**2 - 7729*z**2)


def test_no_square_ternary_quadratic():
    # Commented out test cases should be uncommented after
    # the bug with factor_list() gets merged.

    assert check_solutions(2*x*y + y*z - 3*x*z)
    assert check_solutions(189*x*y - 345*y*z - 12*x*z)
    #assert check_solutions(23*x*y + 34*y*z)
    assert check_solutions(x*y + y*z + z*x)
    assert check_solutions(23*x*y + 23*y*z + 23*x*z)


def test_descent():

    u = ([(13, 23), (3, -11), (41, -113), (91, -3), (1, 1), (1, -1), (17, 13), (123689, 1), (19, -570)])
    for a, b in u:
        w, x, y = descent(a, b)
        assert a*x**2 + b*y**2 == w**2


def test_diophantine():
    # Commented out test cases should be uncommented after
    # the bug with factor_list() gets merged.

    assert check_solutions((x - y)*(y - z)*(z - x))
    assert check_solutions((x - y)*(x**2 + y**2 - z**2))
    assert check_solutions((x - 3*y + 7*z)*(x**2 + y**2 - z**2))
    assert check_solutions((x**2 - 3*y**2 - 1))
    #assert check_solutions(y**2 + 7*x*y)
    #assert check_solutions(x**2 - 3*x*y + y**2)
    #assert check_solutions(z*(x**2 - y**2 - 15))
    #assert check_solutions(x*(2*y - 2*z + 5))
    assert check_solutions((x**2 - 3*y**2 - 1)*(x**2 - y**2 - 15))
    assert check_solutions((x**2 - 3*y**2 - 1)*(y - 7*z))
    assert check_solutions((x**2 + y**2 - z**2)*(x - 7*y - 3*z + 4*w))
    # Following test case caused problems in parametric representation
    # But this can be solved by factroing out y.
    # No need to use methods for ternary quadratic equations.
    #assert check_solutions(y**2 - 7*x*y + 4*y*z)
    assert check_solutions(x**2 - 2*x + 1)


def test_general_pythagorean():

    from sympy.abc import a, b, c, d, e

    assert check_solutions(a**2 + b**2 + c**2 - d**2)
    assert check_solutions(a**2 + 4*b**2 + 4*c**2 - d**2)
    assert check_solutions(9*a**2 + 4*b**2 + 4*c**2 - d**2)
    assert check_solutions(9*a**2 + 4*b**2 - 25*d**2 + 4*c**2 )
    assert check_solutions(9*a**2 - 16*d**2 + 4*b**2 + 4*c**2)
    assert check_solutions(-e**2 + 9*a**2 + 4*b**2 + 4*c**2 + 25*d**2)
    assert check_solutions(16*a**2 - b**2 + 9*c**2 + d**2 + 25*e**2)


def test_diop_general_sum_of_squares():

    from sympy.abc import a, b, c, d, e, f, g, h, i

    assert check_solutions(a**2 + b**2 + c**2 - 5)
    assert check_solutions(a**2 + b**2 + c**2 - 57)
    assert check_solutions(a**2 + b**2 + c**2 - 349560695)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 - 304)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 - 23345)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 - 23345494)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 + e**2 - 1344545)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 + e**2 + f**2 - 6933949)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 + e**2 + f**2 + g**2 - 753934)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 + e**2 + f**2 + g**2 + h**2 - 5)
    assert check_solutions(a**2 + b**2 + c**2 + d**2 + e**2 + f**2 + g**2 + h**2 + i**2 - 693940)


def test_partition():

    tests = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    for test in tests:
        f = partition(test)
        while True:
            try:
                l = next(f)
            except StopIteration:
                break

    tests_k = [8, 10]

    for test in tests_k:
        for k in range(8):
            f = partition(test, k)

            while True:
                try:
                    l = next(f)
                    assert len(l) == k
                except StopIteration:
                    break


def test_prime_as_sum_of_two_squares():

    for i in [5, 13, 17, 29, 37, 41, 2341, 3557, 34841, 64601]:
        a, b = prime_as_sum_of_two_squares(i)
        assert a**2 + b**2 == i


def test_sum_of_three_squares():

    for i in [0, 1, 2, 34, 123, 34304595905, 34304595905394941, 343045959052344,
              800, 801, 802, 803, 804, 805, 806]:
        a, b, c = sum_of_three_squares(i)
        assert a**2 + b**2 + c**2 == i

    assert sum_of_three_squares(7) == (None, None, None)
    assert sum_of_three_squares((4**5)*15) == (None, None, None)


def test_sum_of_four_squares():

    from random import randint

    for i in range(10):
        n = randint(1, 100000000000000)
        a, b, c, d = sum_of_four_squares(n)
        assert a**2 + b**2 + c**2 + d**2 == n


def test_power_representation():

    tests = [(1729, 3, 2), (234, 2, 4), (2, 1, 2), (3, 1, 3), (5, 2, 2), (12352, 2, 4),
             (32760, 2, 3)]

    for test in tests:
        n, p, k = test
        f = power_representation(n, p, k)

        while True:
            try:
                l = next(f)
                assert len(l) == k

                chk_sum = 0
                for l_i in l:
                    chk_sum = chk_sum + l_i**p
                assert chk_sum == n

            except StopIteration:
                break


def test_assumptions():
    """
    Test whether diophantine respects the assumptions.
    """
    #Test case taken from the below so question regarding assumptions in diophantine module
    #http://stackoverflow.com/questions/23301941/how-can-i-declare-natural-symbols-with-sympy
    m, n = symbols('m n', integer=True, positive=True)
    diof = diophantine(n**2 + m * n - 500)
    assert diof == set([(5, 20), (40, 10), (95, 5), (121, 4), (248, 2), (499, 1)])

    a, b = symbols('a b', integer=True, positive=False)
    diof = diophantine(a*b + 2*a + 3*b - 6)
    assert diof == set([(-15, -3), (-9, -4), (-7, -5), (-6, -6), (-5, -8), (-4, -14)])



def check_solutions(eq):
    """
    Determines whether solutions returned by diophantine() satisfy the original
    equation. Hope to generalize this so we can remove functions like check_ternay_quadratic,
    check_solutions_normal, check_solutions()
    """
    s = diophantine(eq)

    terms = factor_list(eq)[1]

    var = list(eq.free_symbols)
    var.sort(key=default_sort_key)

    okay = True

    while len(s) and okay:
        solution = s.pop()

        okay = False

        for term in terms:
            subeq = term[0]

            if simplify(_mexpand(Subs(subeq, var, solution).doit())) == 0:
                okay = True
                break

    return okay

def check_integrality(eq):
    """
    Check that the solutions returned by diophantine() are integers.
    This should be seldom needed except for general quadratic
    equations which are solved with rational transformations.
    """
    def _check_values(x):
        """ Check a number of values. """
        for i in range(-4, 4):
            if not isinstance(simplify(x.subs(t, i)), Integer):
                return False
        return True

    for soln in diophantine(eq, param=t):
        for x in soln:
            if not _check_values(x):
                return False

    return True
