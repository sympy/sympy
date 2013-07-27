from sympy.solvers.diophantine import (diop_solve, diop_pell, diop_bf_pell, length, transformation_to_pell, find_DN, equivalent,
    parametrize_ternary_quadratic, square_factor, pairwise_prime, diop_ternary_quadratic, diop_ternary_quadratic_normal, descent,
    ldescent, classify_diop)

from sympy import symbols, Integer, Matrix, simplify, Subs, S, factorint
from sympy.utilities.pytest import XFAIL, slow

x, y, z, w, t, X, Y = symbols("x, y, z, w, t, X, Y", Integer=True)


def test_linear():

    assert diop_solve(2*x + 3*y - 5) == {x: 3*t - 5, y: -2*t + 5}
    assert diop_solve(3*y + 2*x - 5) == {x: 3*t - 5, y: -2*t + 5}
    assert diop_solve(2*x - 3*y - 5) == {x: -3*t - 5, y: -2*t - 5}
    assert diop_solve(-2*x - 3*y - 5) == {x: -3*t + 5, y: 2*t - 5}
    assert diop_solve(7*x + 5*y) == {x: 5*t, y: -7*t}
    assert diop_solve(2*x + 4*y) == {x: 2*t, y: -t}
    assert diop_solve(4*x + 6*y - 4) == {x: 3*t - 2, y: -2*t + 2}
    assert diop_solve(4*x + 6*y - 3) == {x: None, y: None}
    assert diop_solve(4*x + 3*y -4*z + 5) == \
           {x: 3*t + 4*z - 5, y: -4*t - 4*z + 5, z: z}
    assert diop_solve(4*x + 2*y + 8*z - 5) == {x: None, y: None, z: None}
    assert diop_solve(5*x + 7*y - 2*z - 6) == \
           {x: 7*t + 6*z + 18, y: -5*t - 4*z - 12, z: z}
    assert diop_solve(3*x - 6*y + 12*z - 9) == \
           {x: -2*t - 4*z + 3, y: -t, z: z}
    assert diop_solve(x + 3*y - 4*z + w - 6) == \
           {w: t, x: -t - 3*y + 4*z + 6, y: y, z: z}


def solutions_ok_quadratic(eq):
    """
    Determines whether solutions returned by diop_solve() satisfy the original
    equation.
    """
    s = diop_solve(eq)
    x, y = symbols("x, y", Integer=True)
    ok = True

    while len(s) and ok:
        u, v = s.pop()

        if simplify(simplify(Subs(eq, (x, y), (u, v)).doit())) != 0:
            ok = False
    return ok


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

def test_quadratic_elliptical_case():
    # Elliptical case: B**2 - 4AC < 0
    assert diop_solve(42*x**2 + 8*x*y + 15*y**2 + 23*x + 17*y - 4915) == set([(-Integer(11), -Integer(1))])
    assert diop_solve(4*x**2 + 3*y**2 + 5*x - 11*y + 12) == set([])
    assert diop_solve(x**2 + y**2 + 2*x + 2*y + 2) == set([(-Integer(1), -Integer(1))])
    assert diop_solve(15*x**2 - 9*x*y + 14*y**2 - 23*x - 14*y - 4950) == set([(-Integer(15), Integer(6))])
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

def test_quadratic_non_perfect_square():
    # B**2 - 4*A*C is not a perfect square
    # Used solutions_ok_quadratic() since the solutions are complex expressions involving
    # square roots and exponents
    assert solutions_ok_quadratic(x**2 - 2*x - 5*y**2) == True
    assert solutions_ok_quadratic(3*x**2 - 2*y**2 - 2*x - 2*y) == True
    assert solutions_ok_quadratic(x**2 - x*y - y**2 - 3*y) == True
    assert solutions_ok_quadratic(x**2 - 9*y**2 - 2*x - 6*y) == True

@slow
def test_quadratic_non_perfect_slow():
    assert solutions_ok_quadratic(8*x**2 + 10*x*y - 2*y**2 - 32*x - 13*y - 23) == True
    assert solutions_ok_quadratic(5*x**2 - 13*x*y + y**2 - 4*x - 4*y - 15) == True
    assert solutions_ok_quadratic(-3*x**2 - 2*x*y + 7*y**2 - 5*x - 7) == True

@XFAIL
def test_quadratic_bugs():
    assert diop_solve(x**2 - y**2 - 2*x - 2*y) == set([(t, -t), (-t, -t - 2)])
    assert diop_solve(x**2 - 9*y**2 - 2*x - 6*y) == set([(-3*t + 2, -t), (3*t, -t)])
    assert diop_solve(4*x**2 - 9*y**2 - 4*x - 12*y - 3) == set([(-3*t - 3, -2*t - 3), (3*t + 1, -2*t - 1)])


def test_pell():
    # Most of the test cases were adapted from,
    # Solving the generalized Pell equation x**2 - D*y**2 = N, John P. Robertson, July 31, 2004.
    # http://www.jpr2718.org/pell.pdf
    # others are verified using Wolfram Alpha.

    # Covers cases where D <= 0 or D > 0 and D is a square or N = 0
    # Solutions are straightforward in these cases.
    assert diop_pell(3, 0) == [(0, 0)]
    assert diop_pell(-17, -5) == []
    assert diop_pell(-19, 23) == [(2, 1)]
    assert diop_pell(-13, 17) == [(2, 1)]
    assert diop_pell(-15, 13) == []
    assert diop_pell(0, 5) == []
    assert diop_pell(0, 9) == [(3, t), (-3, t)]
    assert diop_pell(9, 0) == [(3*t, t), (-3*t, t)]
    assert diop_pell(16, 24) == []
    assert diop_pell(9, 180) == [(18, 4)]
    assert diop_pell(9, -180) == [(12, 6)]
    assert diop_pell(7, 0) == [(0, 0)]

    # D > 0 and D is not a square

    # N = 1
    assert diop_pell(13, 1) == [(649, 180)]
    assert diop_pell(980, 1) == [(51841, 1656)]
    assert diop_pell(981, 1) == [(158070671986249, 5046808151700)]
    assert diop_pell(986, 1) == [(49299, 1570)]
    assert diop_pell(991, 1) == [(379516400906811930638014896080, 12055735790331359447442538767)]
    assert diop_pell(17, 1) == [(33, 8)]
    assert diop_pell(19, 1) == [(170, 39)]

    # N = -1
    assert diop_pell(13, -1) == [(18, 5)]
    assert diop_pell(991, -1) == []
    assert diop_pell(41, -1) == [(32, 5)]
    assert diop_pell(290, -1) == [(17, 1)]
    assert diop_pell(21257, -1) == [(13913102721304, 95427381109)]
    assert diop_pell(32, -1) == []

    # |N| > 1
    # Some tests were created using calculator at
    # http://www.numbertheory.org/php/patz.html

    assert diop_pell(13, -4) == [(3, 1), (393, 109), (36, 10)]
    # Source I referred returned (3, 1), (393, 109) and (-3, 1) as fundamental solutions
    # So (-3, 1) and (393, 109) should be in the same equivalent class
    assert equivalent(-3, 1, 393, 109, 13, -4) == True

    assert diop_pell(13, 27) == [(220, 61), (40, 11), (768, 213), (12, 3)]
    assert set(diop_pell(157, 12)) == \
    set([(Integer(13), Integer(1)), (Integer(10663), Integer(851)), (Integer(579160), Integer(46222)), \
        (Integer(483790960),Integer(38610722)), (Integer(26277068347), Integer(2097138361)), (Integer(21950079635497), Integer(1751807067011))])
    assert diop_pell(13, 25) == [(3245, 900)]
    assert diop_pell(192, 18) == []
    assert diop_pell(23, 13) == [(-6, 1), (6, 1)]
    assert diop_pell(167, 2) == [(13, 1)]
    assert diop_pell(167, -2) == []

    assert diop_pell(123, -2) == [(11, 1)]
    # One calculator returned [(11, 1), (-11, 1)] but both of these are in
    # the same equivalence class
    assert equivalent(11, 1, -11, 1, 123, -2) == True

    assert diop_pell(123, -23) == [(-10, 1), (10, 1)]


def test_bf_pell():

    assert diop_bf_pell(13, -4) == [(3, 1), (-3, 1), (36, 10)]
    assert diop_bf_pell(13, 27) == [(12, 3), (-12, 3), (40, 11), (-40, 11)]
    assert diop_bf_pell(167, -2) == []
    assert diop_bf_pell(1729, 1) == [(44611924489705, 1072885712316)]
    assert diop_bf_pell(89, -8) == [(9, 1), (-9, 1)]
    assert diop_bf_pell(21257, -1) == [(13913102721304, 95427381109)]
    assert diop_bf_pell(340, -4) == [(756, 41)]


def test_length():

    assert length(-2, 4, 5) == 3
    assert length(-5, 4, 17) == 4
    assert length(0, 4, 13) == 6
    assert length(-31, 8, 613) == 67
    assert length(7, 13, 11) == 23
    assert length(-40, 5, 23) == 4


def is_transformation_ok(eq):
    """
    Test whether X*Y, X, or Y terms are present in the equation
    after transforming the equation using the transformation returned
    by transformation_to_pell(). If they are not present we are good.
    Moreover, coefficient of X**2 should be a divisor of coefficient of
    Y**2 and the constant term.
    """
    A, B = transformation_to_pell(eq)
    u = (A*Matrix([X, Y]) + B)[0]
    v = (A*Matrix([X, Y]) + B)[1]
    simplified = simplify(Subs(eq, (x, y), (u, v)).doit())

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

    assert is_transformation_ok(-13*x**2 - 7*x*y + y**2 + 2*x - 2*y - 14) == True
    assert is_transformation_ok(-17*x**2 + 19*x*y - 7*y**2 - 5*x - 13*y - 23) == True
    assert is_transformation_ok(x**2 - y**2 + 17) == True
    assert is_transformation_ok(-x**2 + 7*y**2 - 23) == True
    assert is_transformation_ok(25*x**2 - 45*x*y + 5*y**2 - 5*x - 10*y + 5) == True
    assert is_transformation_ok(190*x**2 + 30*x*y + y**2 - 3*y - 170*x - 130) == True
    assert is_transformation_ok(x**2 - 2*x*y -190*y**2 - 7*y - 23*x - 89) == True


def test_find_DN():

    # Note that b**2 - 4*a*c > 0 and should not be a perfect square for this
    # method to work
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


def check_ternary_quadratic_normal(eq):
    var, coeff, type = classify_diop(eq)
    x_0, y_0, z_0 = diop_ternary_quadratic_normal(eq)

    x = var[0]
    y = var[1]
    z = var[2]

    return (x_0**2*coeff[x**2] + y_0**2*coeff[y**2] + z_0**2*coeff[z**2] == 0)


def test_diop_ternary_quadratic_normal():

    assert check_ternary_quadratic_normal(234*x**2 - 65601*y**2 - z**2) == True
    assert check_ternary_quadratic_normal(23*x**2 + 616*y**2 - z**2) == True
    assert check_ternary_quadratic_normal(5*x**2 + 4*y**2 - z**2) == True
    assert check_ternary_quadratic_normal(3*x**2 + 6*y**2 - 3*z**2) == True
    assert check_ternary_quadratic_normal(x**2 + 3*y**2 - z**2) == True
    assert check_ternary_quadratic_normal(4*x**2 + 5*y**2 - z**2) == True
    assert check_ternary_quadratic_normal(x**2 + y**2 - z**2) == True
    assert check_ternary_quadratic_normal(16*x**2 + y**2 - 25*z**2) == True
    assert check_ternary_quadratic_normal(6*x**2 - y**2 + 10*z**2) == True
    assert check_ternary_quadratic_normal(213*x**2 + 12*y**2 - 9*z**2) == True
    assert check_ternary_quadratic_normal(34*x**2 - 3*y**2 - 301*z**2) == True
    assert check_ternary_quadratic_normal(124*x**2 - 30*y**2 - 7729*z**2) == True


def check_ternary_quadratic(eq):
    var, coeff, diop_type = classify_diop(eq)
    x_0, y_0, z_0 = diop_ternary_quadratic(eq)

    x = var[0]
    y = var[1]
    z = var[2]

    return (x_0**2*coeff[x**2] + y_0**2*coeff[y**2] + z_0**2*coeff[z**2] + x_0*y_0*coeff[x*y]
            + y_0*z_0*coeff[y*z] + z_0*x_0*coeff[z*x] == 0)

def test_diop_ternary_quadratic():

    assert check_ternary_quadratic(2*x**2 + z**2 + y**2 - 4*x*y) == True
    assert check_ternary_quadratic(x**2 - y**2 - z**2 - x*y - y*z) == True
    assert check_ternary_quadratic(3*x**2 - x*y - y*z - x*z) == True
    assert check_ternary_quadratic(x**2 - y*z - x*z) == True
    assert check_ternary_quadratic(5*x**2 - 3*x*y - x*z) == True
    assert check_ternary_quadratic(4*x**2 - 5*y**2 - x*z) == True
    assert check_ternary_quadratic(3*x**2 + 2*y**2 - z**2 - 2*x*y + 5*y*z - 7*y*z) == True
    assert check_ternary_quadratic(8*x**2 - 12*y*z) == True
    assert check_ternary_quadratic(45*x**2 - 7*y**2 - 8*x*y - z**2) == True
    assert check_ternary_quadratic(x**2 - 49*y**2 - z**2 + 13*z*y -8*x*y) == True
    assert check_ternary_quadratic(90*x**2 + 3*y**2 + 5*x*y + 2*z*y + 5*x*z) == True

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


def check_parametrize_ternary_quadratic(eq):

    x_p, y_p, z_p = parametrize_ternary_quadratic(eq)
    var, jnk1, jnk2 = classify_diop(eq)

    x = var[0]
    y = var[1]
    z = var[2]

    return simplify(simplify(Subs(eq, (x, y, z), (x_p, y_p, z_p)).doit())) == 0


def test_parametrize_ternary_quadratic():

    assert check_parametrize_ternary_quadratic(x**2 + y**2 - z**2) == True
    assert check_parametrize_ternary_quadratic(x**2 + 2*x*y + z**2) == True
    assert check_parametrize_ternary_quadratic(234*x**2 - 65601*y**2 - z**2) == True
    assert check_parametrize_ternary_quadratic(3*x**2 + 2*y**2 - z**2 - 2*x*y + 5*y*z - 7*y*z) == True
    assert check_parametrize_ternary_quadratic(x**2 - y**2 - z**2) == True
    assert check_parametrize_ternary_quadratic(x**2 - 49*y**2 - z**2 + 13*z*y - 8*x*y) == True
    assert check_parametrize_ternary_quadratic(8*x*y + z**2) == True
    assert check_parametrize_ternary_quadratic(124*x**2 - 30*y**2 - 7729*z**2) == True
    assert check_parametrize_ternary_quadratic(236*x**2 - 225*y**2 - 11*x*y - 13*y*z - 17*x*z) == True
    assert check_parametrize_ternary_quadratic(90*x**2 + 3*y**2 + 5*x*y + 2*z*y + 5*x*z) == True


def test_descent():

    u = ([(13, 23), (3, -11), (41, -113), (91, -3), (1, 1), (1, -1), (17, 13), (123689, 1), (19, -570)])
    for a, b in u:
        w, x, y = descent(a, b)
        assert a*x**2 + b*y**2 == w**2


def test_bug_parametrize_ternary_quadratic():

    parametrize_ternary_quadratic(y**2 - 7*x*y + 4*y*z)
