from sympy.core.symbol import symbols
from sympy.simplify.simplify import simplify
from sympy.core.function import expand
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.series.taylor_2D import TaylorTwovariable

X, Y = symbols("X Y")

def test_poly_simple():
    X, Y = symbols("X Y")
    result = TaylorTwovariable("X**2 + Y**2", 0, 0, 2)
    expected = X**2 + Y**2
    assert simplify(result - expected) == 0

def test_poly_mixed():
    X, Y = symbols("X Y")
    result = TaylorTwovariable("X*Y + X**2", 1, 1, 2)
    expected = 2 + 3*(X-1) + 1*(Y-1) + (X-1)**2 + (X-1)*(Y-1)
    assert simplify(result - expected) == 0

def test_sin():
    X, Y = symbols("X Y")
    result = TaylorTwovariable("sin(X+Y)", 0, 0, 3)
    expected = (X + Y) - (X + Y)**3/6
    assert simplify(result - expected) == 0

def test_cos_case():
    result = simplify(TaylorTwovariable("cos(X*Y)", 0, 0, 4))
    assert simplify(result - (1 - (X*Y)**2/2)) == 0

def test_exponential_case():
    result = simplify(TaylorTwovariable("exp(X + Y)", 0, 0, 3))
    expected = expand(1 + (X+Y) + (X+Y)**2/2 + (X+Y)**3/6)
    assert simplify(result - expected) == 0

def test_constant_case():
    assert simplify(TaylorTwovariable("5", 0, 0, 3) - 5) == 0

def test_offset_case():
    result = simplify(TaylorTwovariable("X**2 + Y**2", 1, 1, 3))
    expected = expand(2 + 2*(X-1) + 2*(Y-1) + (X-1)**2 + (Y-1)**2)
    assert simplify(result - expected) == 0


def test_log_case():
    result = TaylorTwovariable("log(1 + X + Y)", 0, 0, 3)
    approx = expand((X+Y) - (X+Y)**2/2 + (X+Y)**3/3)
    assert simplify(result - approx) == 0

def test_asymmetric_case():
    result = TaylorTwovariable("X**3 + 2*Y", 0, 0, 3)
    assert simplify(result - (X**3 + 2*Y)) == 0

def test_higher_order_case():
    result = TaylorTwovariable("exp(X*Y)", 0, 0, 4)
    expected = expand((X**2*Y**2)/2 + X*Y + 1)
    assert simplify(result - expected) == 0

def test_shifted_trig():
    result = TaylorTwovariable("sin(X+Y)", 1, 1, 2)
    approx = expand((sin(2) + cos(2)*(X-1) + cos(2)*(Y-1) - sin(2)*((X-1)+(Y-1))**2/2))
    assert simplify(result - approx) == 0
