from sympy.solvers.decompogen import decompogen
from sympy import sin, cos, sqrt, Abs
from sympy import symbols
from sympy.utilities.pytest import XFAIL

x, y = symbols('x y')

def test_decompogen():
    assert decompogen(sin(cos(x)), x) == [sin(x), cos(x)]
    assert decompogen(sin(x)**2 + sin(x) + 1, x) == [x**2 + x + 1, sin(x)]
    assert decompogen(sqrt(6*x**2 - 5), x) == [sqrt(x), 6*x**2 - 5]
    assert decompogen(sin(sqrt(cos(x**2 + 1))), x) == [sin(x), sqrt(x), cos(x), x**2 + 1]
    assert decompogen(Abs(cos(x)**2 + 3*cos(x) - 4), x) == [Abs(x), x**2 + 3*x - 4, cos(x)]
    assert decompogen(sin(x)**2 + sin(x) - sqrt(3)/2, x) == [x**2 + x - sqrt(3)/2, sin(x)]
    assert decompogen(Abs(cos(y)**2 + 3*cos(x) - 4), x) == [Abs(x), 3*x + cos(y)**2 - 4, cos(x)]

def test_decompogen_poly():
    assert decompogen(x**4 + 2*x**2 + 1, x) == [x**2 + 2*x + 1, x**2]
    assert decompogen(x**4 + 2*x**3 - x - 1, x) == [x**2 - x - 1, x**2 + x]


@XFAIL
def test_decompogen_fails():
    A = lambda x: x**2 + 2*x + 3
    B = lambda x: 4*x**2 + 5*x + 6
    assert decompogen(A(x*exp(x)), x) == [x**2 + 2*x + 3, x*exp(x)]
    assert decompogen(A(B(x)), x) == [x**2 + 2*x + 3, 4*x**2 + 5*x + 6]
    assert decompogen(A(1/x + 1/x**2), x) == [x**2 + 2*x + 3, 1/x + 1/x**2]
    assert decompogen(A(1/x + 2/(x + 1)), x) == [x**2 + 2*x + 3, 1/x + 2/(x + 1)]
