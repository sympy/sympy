from sympy import symbols, sin, Matrix
from sympy.printing.lambdarepr import lambdarepr
x,y,z = symbols("xyz")

def test_basic():
    assert lambdarepr(x*y)=="x*y"
    assert lambdarepr(x+y) in ["y + x", "x + y"]
    assert lambdarepr(x**y)=="x**y"

def test_matrix():
    A = Matrix([[x,y],[y*x,z**2]])
    assert lambdarepr(A)=="Matrix([[  x,    y],[x*y, z**2]])"