########################################### TESTS ##############################################################################
from sympy import symbols
from sympy import Integer
from sympy.solvers.diophantine import diop_solve, diop_pell, diop_bf_pell, length
x, y, z, w, t = symbols("x, y, z, w, t", integer=True)


def test_linear():

    assert diop_solve(2*x + 3*y - 5) == {x: 15*t - 5, y: -10*t + 5}
    assert diop_solve(3*y + 2*x - 5) == {x: 15*t - 5, y: -10*t + 5}
    assert diop_solve(2*x - 3*y - 5) == {x: -15*t - 5, y: -10*t - 5}
    assert diop_solve(-2*x - 3*y - 5) == {x: -15*t + 5, y: 10*t - 5}
    assert diop_solve(7*x + 5*y) == {x: 5*t, y: -7*t}
    assert diop_solve(2*x + 4*y) == {x: 2*t, y: -t}
    assert diop_solve(4*x + 6*y - 4) == {x: 6*t - 2, y: -4*t + 2}
    assert diop_solve(4*x + 6*y - 3) == {x: None, y: None}
    assert diop_solve(4*x + 3*y -4*z + 5) == \
           {x: -15*t + 4*z - 5, y: 20*t - 4*z + 5, z: z}
    assert diop_solve(4*x + 2*y + 8*z - 5) == {x: None, y: None, z: None}
    assert diop_solve(5*x + 7*y - 2*z - 6) == \
           {x: 42*t + 6*z + 18, y: -30*t - 4*z - 12, z: z}
    assert diop_solve(3*x - 6*y + 12*z - 9) == \
           {x: -6*t - 4*z + 3, y: -3*t, z: z}
    assert diop_solve(x + 3*y - 4*z + w - 6) == \
           {w: 6*t, x: -6*t - 3*y + 4*z + 6, y: y, z: z}


def test_quadratic():

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

    # Elliptical case: B**2 - 4AC < 0
    assert diop_solve(42*x**2 + 8*x*y + 15*y**2 + 23*x + 17*y - 4915) == set([(-Integer(11), -Integer(1))])
    assert diop_solve(4*x**2 + 3*y**2 + 5*x - 11*y + 12) == set([])
    assert diop_solve(x**2 + y**2 + 2*x + 2*y + 2) == set([(-Integer(1), -Integer(1))])
    assert diop_solve(15*x**2 - 9*x*y + 14*y**2 - 23*x - 14*y - 4950) == set([(-Integer(15), Integer(6))])
    assert diop_solve(10*x**2 + 12*x*y + 12*y**2 - 34) == set([(Integer(1), -Integer(2)), (-Integer(1), -Integer(1)),\
        (Integer(1), Integer(1)), (-Integer(1), Integer(2))])
    assert diop_solve(3*x**2 + 5*x*y + 7*y**2) == set([(Integer(0), Integer(0))])

    # Parabolic case: B**2 - 4AC = 0
    assert diop_solve(8*x**2 - 24*x*y + 18*y**2 + 5*x + 7*y + 16) == \
        set([(-174*t**2 + 17*t - 2, -116*t**2 + 21*t - 2), (-174*t**2 + 41*t - 4, -116*t**2 + 37*t - 4)])
    assert diop_solve(8*x**2 - 24*x*y + 18*y**2 + 6*x + 12*y - 6) == \
        set([(-63*t**2 + 12*t, -42*t**2 + 15*t -1), (-63*t**2 + 30*t - 3, -42*t**2 + 27*t - 4)])
    assert diop_solve(8*x**2 + 24*x*y + 18*y**2 + 4*x + 6*y - 7) == \
        set([])
    assert diop_solve(x**2 + 2*x*y + y**2 + 2*x + 2*y + 1) == set([(-t,t - 1)])
    assert diop_solve(x**2 - 2*x*y + y**2 + 2*x + 2*y + 1) == \
        set([(-4*t**2, -4*t**2 + 4*t - 1),(-4*t**2 + 4*t -1, -4*t**2 + 8*t - 4)])

    # Random
    assert diop_solve(48*x*y) == set([(Integer(0), t), (t, Integer(0))])
    assert diop_solve(4*x**2 - 5*x*y + y**2 + 2) == \
        set([(-Integer(1), -Integer(3)),(-Integer(1), -Integer(2)),(Integer(1), Integer(2)),(Integer(1), Integer(3))])


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

    # D > 0 and D is square free

    # N == 1
    assert diop_pell(13, 1) == [(649, 180)]
    assert diop_pell(980, 1) == [(51841, 1656)]
    assert diop_pell(981, 1) == [(158070671986249, 5046808151700)]
    assert diop_pell(986, 1) == [(49299, 1570)]
    assert diop_pell(991, 1) == [(379516400906811930638014896080, 12055735790331359447442538767)]
    assert diop_pell(17, 1) == [(33, 8)]

    # N == -1
    assert diop_pell(13, -1) == [(18, 5)]
    assert diop_pell(991, -1) == []
    assert diop_pell(41, -1) == [(32, 5)]
    assert diop_pell(290, -1) == [(17, 1)]
    assert diop_pell(21257, -1) == [(13913102721304, 95427381109)]
    assert diop_pell(32, -1) == []

    # |N| > 1, need more tests
    assert diop_pell(13, -4) == [(3, 1), (393, 109), (36, 10)] # check this
    assert diop_pell(13, 27) == [(220, 61), (40, 11), (768, 213), (12, 3)]
    assert set(diop_pell(157, 12)) == \
    set([(Integer(13), Integer(1)), (Integer(10663), Integer(851)), (Integer(579160), Integer(46222)), \
        (Integer(483790960),Integer(38610722)), (Integer(26277068347), Integer(2097138361)), (Integer(21950079635497), Integer(1751807067011))])


def test_diop_bf_pell():

    assert diop_bf_pell(13, -4) == [(3, 1), (-3, 1), (36, 10)]
    assert diop_bf_pell(13, 27) == [(12, 3), (-12, 3), (40, 11), (-40, 11)]
    # More tests should be added


def test_length():

    assert length(-2, 4, 5) == 3
    assert length(-5, 4, 17) == 4
    assert length(0, 4, 13) == 6
    assert length(-31, 8, 613) == 67
    assert length(7, 13, 11) == 23
    assert length(-40, 5, 23) == 4
