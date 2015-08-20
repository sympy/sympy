from sympy.solvers.decompogen import decompogen
from sympy import sin, cos, sqrt, Abs
from sympy import Symbol

x = Symbol('x')


def test_decompogen():
    assert decompogen(sin(cos(x)), x) == [sin(x), cos(x)]
    assert decompogen(sin(x)**2 + sin(x) + 1, x) == [x**2 + x + 1, sin(x)]
    assert decompogen(sqrt(6*x**2 - 5), x) == [sqrt(x), 6*x**2 - 5]
    assert decompogen(sin(sqrt(cos(x**2 + 1))), x) == [sin(x), sqrt(x), cos(x), x**2 + 1]
    assert decompogen(Abs(cos(x)**2 + 3*cos(x) - 4), x) == [Abs(x), x**2 + 3*x - 4, cos(x)]
