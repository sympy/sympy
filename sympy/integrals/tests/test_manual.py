from sympy import (sin, cos, tan, sec, csc, cot, log, exp, Abs, atan,
                   Symbol)
from sympy.integrals.manualintegrate import manualintegrate, integral_steps

x = Symbol('x')
y = Symbol('y')

def test_manualintegrate_polynomials():
    assert manualintegrate(y, x) == x*y
    assert manualintegrate(exp(2), x) == x * exp(2)
    assert manualintegrate(x**2, x) == x**3 / 3
    assert manualintegrate(3 * x**2 + 4 * x**3, x) == x**3 + x**4

    assert manualintegrate((x + 2)**3, x) == (x + 2)**4 / 4
    assert manualintegrate((3*x + 4)**2, x) == (3*x + 4)**3 / 9

def test_manualintegrate_exponentials():
    assert manualintegrate(exp(2*x), x) == exp(2*x) / 2
    assert manualintegrate(2**x, x) == (2 ** x) / log(2)

    assert manualintegrate(1 / x, x) == log(Abs(x))
    assert manualintegrate(1 / (2*x + 3), x) == log(Abs(2*x + 3)) / 2

def test_manualintegrate_trigonometry():
    assert manualintegrate(sin(x), x) == -cos(x)
    assert manualintegrate(tan(x), x) == -log(Abs(cos(x)))

#    assert manualintegrate(sec(x), x) == log(Abs(sec(x) + tan(x)))
#    assert manualintegrate(csc(x), x) == -log(Abs(csc(x) + cot(x)))

    assert manualintegrate(sin(x) * cos(x), x) == -cos(x) ** 2 / 2

def test_manualintegrate_trigpowers():
    assert manualintegrate(sin(x)**2 * cos(x), x) == sin(x)**3 / 3
    assert manualintegrate(sin(x)**2 * cos(x) **2, x) == \
        x / 8 - sin(4*x) / 32
    assert manualintegrate(sin(x) * cos(x)**3, x) == -cos(x)**4 / 4
    assert manualintegrate(sin(x)**3 * cos(x)**2, x) == \
        -cos(x)**5 / 5 + cos(x)**3 / 3

    assert manualintegrate(tan(x)**3 * sec(x), x) == sec(x)**3/3 - sec(x)
    assert manualintegrate(tan(x) * sec(x) **2, x) == sec(x)**2/2

    # print integral_steps(cot(x) ** 5 * csc(x), x)
    # print integral_steps(cot(x)**2 * csc(x)**6, x)

def test_manualintegrate_inversetrig():
    assert manualintegrate(exp(x) / (1 + exp(2*x)), x) == atan(exp(x))
    assert manualintegrate(1 / (4 + 9 * x**2), x) == atan(3 * x/2) / 6
    assert manualintegrate(1 / (16 + 16 * x**2), x) == atan(x) / 16
    assert manualintegrate(1 / (4 + x**2), x) == atan(x / 2) / 2
    assert manualintegrate(1 / (1 + 4 * x**2), x) == atan(2*x) / 2
