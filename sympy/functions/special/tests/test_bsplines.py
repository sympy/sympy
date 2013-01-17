from sympy.functions import bspline_basis, bspline_basis_set
from sympy import Piecewise, Interval
from sympy import symbols, Rational

x, y = symbols('x,y')


def test_basic_degree_0():
    d = 0
    knots = range(5)
    splines = bspline_basis_set(d, knots, x)
    for i in range(len(splines)):
        assert splines[i] == Piecewise((1, Interval(i, i + 1)
                                       .contains(x)), (0, True))


def test_basic_degree_1():
    d = 1
    knots = range(5)
    splines = bspline_basis_set(d, knots, x)
    assert splines[0] == Piecewise(
        (x, Interval(0, 1, False, True).contains(x)),
        (2 - x, Interval(1, 2).contains(x)), (0, True))
    assert splines[1] == Piecewise(
        (-1 + x, Interval(1, 2, False, True).contains(x)),
        (3 - x, Interval(2, 3).contains(x)), (0, True))
    assert splines[2] == Piecewise(
        (-2 + x, Interval(2, 3, False, True).contains(x)),
        (4 - x, Interval(3, 4).contains(x)), (0, True))


def test_basic_degree_2():
    d = 2
    knots = range(5)
    splines = bspline_basis_set(d, knots, x)
    b0 = Piecewise((x**2/2, Interval(0, 1, False, True).contains(x)),
        (Rational(
            -3, 2) + 3*x - x**2, Interval(1, 2, False, True).contains(x)),
        (Rational(9, 2) - 3*x + x**2/2, Interval(2, 3).contains(x)), (0, True))
    b1 = Piecewise(
        (Rational(1, 2) - x + x**2/2, Interval(1, 2, False, True).contains(x)),
        (Rational(
            -11, 2) + 5*x - x**2, Interval(2, 3, False, True).contains(x)),
        (8 - 4*x + x**2/2, Interval(3, 4).contains(x)), (0, True))
    assert splines[0] == b0
    assert splines[1] == b1


def test_basic_degree_3():
    d = 3
    knots = range(5)
    splines = bspline_basis_set(d, knots, x)
    b0 = Piecewise(
        (x**3/6, Interval(0, 1, False, True).contains(x)),
        (Rational(2, 3) - 2*x + 2*x**2 - x**3/2, Interval(1, 2,
         False, True).contains(x)),
        (Rational(-22, 3) + 10*x - 4*x**2 + x**3/2, Interval(2, 3,
         False, True).contains(x)),
        (Rational(32, 3) - 8*x + 2*x**2 - x**3/6, Interval(3, 4).contains(x)),
        (0, True)
    )
    assert splines[0] == b0


def test_repeated_degree_1():
    d = 1
    knots = [0, 0, 1, 2, 2, 3, 4, 4]
    splines = bspline_basis_set(d, knots, x)
    assert splines[0] == Piecewise((1 - x, Interval(0, 1).contains(x)),
                                   (0, True))
    assert splines[1] == Piecewise(
        (x, Interval(0, 1, False, True).contains(x)),
        (2 - x, Interval(1, 2).contains(x)), (0, True))
    assert splines[2] == Piecewise((-1 + x, Interval(1, 2).contains(x)
                                   ), (0, True))
    assert splines[3] == Piecewise((3 - x, Interval(2, 3).contains(x)),
                                   (0, True))
    assert splines[4] == Piecewise(
        (-2 + x, Interval(2, 3, False, True).contains(x)),
        (4 - x, Interval(3, 4).contains(x)), (0, True))
    assert splines[5] == Piecewise((-3 + x, Interval(3, 4).contains(x)
                                   ), (0, True))
