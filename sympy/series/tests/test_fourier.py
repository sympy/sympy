from sympy import (FourierSeries, symbols, pi, Piecewise, sin, cos, Rational,
                   oo, fourier_series)
from sympy.utilities.pytest import raises

x, y, z = symbols('x y z')

fo = FourierSeries(x, (x, -pi, pi))
fe = FourierSeries(x**2, (x, -pi, pi))
fp = FourierSeries(Piecewise((0, x < 0), (pi, True)), (x, -pi, pi))


def test_FourierSeries():
    assert fo.function == x
    assert fo.x == x
    assert fo.period == (-pi, pi)

    assert fo.term(3) == 2*sin(3*x) / 3
    assert fe.term(3) == -4*cos(3*x) / 9
    assert fp.term(3) == 2*sin(3*x) / 3

    assert fo.as_leading_term(x) == 2*sin(x)
    assert fe.as_leading_term(x) == pi**2 / 3
    assert fp.as_leading_term(x) == pi / 2

    assert fo.as_series() == 2*sin(x) - sin(2*x) + (2*sin(3*x) / 3)
    assert fe.as_series() == -4*cos(x) + cos(2*x) + pi**2 / 3
    assert fp.as_series() == 2*sin(x) + (2*sin(3*x) / 3) + pi / 2

    def _check_iter(f, i):
        for ind, t in enumerate(f):
            assert t == f[ind]
            if ind == i:
                break

    _check_iter(fo, 3)
    _check_iter(fe, 3)
    _check_iter(fp, 3)

    raises(ValueError, lambda: FourierSeries(x, (0, 1, 2)))
    raises(ValueError, lambda: FourierSeries(x, (x, 0, oo)))


def test_FourierSeries_2():
    p = Piecewise((0, x < 0), (x, True))
    f = FourierSeries(p, (x, -2, 2))

    assert f.term(3) == (2*sin(3*pi*x / 2) / (3*pi) -
                         4*cos(3*pi*x / 2) / (9*pi**2))
    assert f.as_series() == (2*sin(pi*x / 2) / pi - sin(pi*x) / pi
                      - 4*cos(pi*x / 2) / pi**2 + Rational(1, 2))


def test_fourier_series():
    assert isinstance(fourier_series(x), FourierSeries)


def test_FourierSeries__operations():
    fes = fe.scale(-1).shift(pi**2)
    assert fes.as_series() == 4*cos(x) - cos(2*x) + 2*pi**2 / 3

    assert fp.shift(-pi/2).as_series() == (2*sin(x) + (2*sin(3*x) / 3) +
                                           (2*sin(5*x) / 5))


    fos = fo.scale(3)
    assert fos.as_series() == 6*sin(x) - 3*sin(2*x) + 2*sin(3*x)

    fx = fe.scalex(2).shiftx(1)
    assert fx.as_series() == -4*cos(2*x + 2) + cos(4*x + 4) + pi**2 / 3

    fl = fe.scalex(3).shift(-pi).scalex(2).shiftx(1).scale(4)
    assert fl.as_series() == (-16*cos(6*x + 6) + 4*cos(12*x + 12)
                              + -4*pi + 4*pi**2 / 3)

    raises(ValueError, lambda: fo.shift(x))
    raises(ValueError, lambda: fo.shiftx(sin(x)))
    raises(ValueError, lambda: fo.scale(x*y))
    raises(ValueError, lambda: fo.scalex(x**2))


def test_FourierSeries__neg():
    assert (-fo).as_series() == -2*sin(x) + sin(2*x) - (2*sin(3*x) / 3)
    assert (-fe).as_series() == +4*cos(x) - cos(2*x) - pi**2 / 3


def test_FourierSeries__add__sub():
    assert fo + fo == fo.scale(2)
    assert fo - fo == 0
    assert -fe - fe == fe.scale(-2)

    assert (fo + fe).as_series() == 2*sin(x) - sin(2*x) - 4*cos(x) + cos(2*x) \
        + pi**2 / 3
    assert (fo - fe).as_series() == 2*sin(x) - sin(2*x) + 4*cos(x) - cos(2*x) \
        - pi**2 / 3
