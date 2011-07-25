from sympy import symbols, sin, Matrix, Interval, Piecewise
from sympy.utilities.pytest import raises

from sympy.printing.lambdarepr import lambdarepr

x,y,z = symbols("x,y,z")

def test_basic():
    assert lambdarepr(x*y)=="x*y"
    assert lambdarepr(x+y) in ["y + x", "x + y"]
    assert lambdarepr(x**y)=="x**y"

def test_matrix():
    A = Matrix([[x,y],[y*x,z**2]])
    assert lambdarepr(A)=="Matrix([[  x,    y],[x*y, z**2]])"

def test_piecewise():
    p = Piecewise(
        (x, x < 1),
        (x**2, Interval(3, 4, True, False)),
        (0, True)
    )
    assert lambdarepr(p) == 'iff(x < 1,x,iff(((x <= 4) and (3 < x)),x**2,iff(True,0,0)))'

def test_settings():
    raises(TypeError, 'lambdarepr(sin(x),method="garbage")')
