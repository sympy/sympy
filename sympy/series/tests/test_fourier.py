from sympy import pi, sin, cos, Piecewise, Rational
from sympy.series.fourier import SeqFormula, FourierSeries, fourier_series
from sympy.abc import x, n
from sympy.utilities.pytest import raises
from sympy.core.compatibility import range

def test_SeqFormula_generator():
    seq = SeqFormula(x*(x+1), x, (0, 5))
    assert [val for val in seq] == [n*(n+1) for n in range(6)]

    seq2 = SeqFormula('x*(x+1)', 'x', (0, 5))
    assert [val for val in seq2] == [n*(n+1) for n in range(6)]

def test_SeqFormula_coeff():
    seq = SeqFormula(x*(x+1), x, (0, 5))
    assert seq.coeff(3) == 12

    raises(IndexError, lambda: seq.coeff(10))
    raises(IndexError, lambda: seq.coeff(-1))

def test_FourierSeries_even_function():
    f = FourierSeries(x**2, x)
    assert f.as_nseries(n=3) == -4*cos(x) + cos(2*x) - 4*cos(3*x)/9 \
            + pi**2/3
    assert f[0] == pi**2/3
    assert f[1] == -4*cos(x)

def test_FourierSeries_odd_function():
    f = FourierSeries(x, x)
    assert f.as_nseries(n=3) == 2*sin(x) - sin(2*x) + 2*sin(3*x)/3
    assert f[0] == 0
    assert f[1] == 2*sin(x)

def test_FourierSeries_simple_functions():
    pe = Piecewise((0, x <= 0), (x, True))
    f = FourierSeries(pe, x, (-2, 2))
    assert f.as_nseries(n=3) == 2*sin(pi*x/2)/pi - sin(pi*x)/pi \
            + 2*sin(3*pi*x/2)/(3*pi) - 4*cos(pi*x/2)/pi**2 \
            - 4*cos(3*pi*x/2)/(9*pi**2) + Rational(1,2)
    assert f[0] == Rational(1,2)
    assert f[1] == 2*sin(pi*x/2)/pi - 4*cos(pi*x/2)/pi**2

def test_fourier_series_nterms():
    f = FourierSeries(x, x)
    fn = fourier_series(x, x, n=3)
    assert f.as_nseries(n=3) == fn

def test_fourier_series_generator():
    f = FourierSeries(x, x)
    fn = fourier_series(x, x, n=None)
    assert next(fn) == f[0]
    assert next(fn) == f[1]
    assert next(fn) == f[2]
