from sympy.series.kauers import finite_diff
from sympy.abc import x, y, z, w, n
from sympy import sin, cos
from sympy import pi
from sympy.series.kauers import dominant_term

def test_finite_diff():
    assert finite_diff(x**2 + 2*x + 1, x) == 2*x + 3
    assert finite_diff(y**3 + 2*y**2 + 3*y +5, y) == 3*y**2 + 7*y + 6
    assert finite_diff(z**2 - 2*z + 3, z) == 2*z - 1
    assert finite_diff(w**2 + 3*w -2, w) == 2*w + 4
    assert finite_diff(sin(x), x,  pi/6) == -sin(x) + sin(x + pi/6)
    assert finite_diff(cos(y),y,  pi/3)
    assert finite_diff(x**2 - 2*x + 3, x,  2) == 4*x
    assert finite_diff(n**2 - 2*n + 3, n, 3) == 6*n + 3

def test_dominant_term():
    assert dominant_term(3*n**2 + 4*n + 8, n) ==  n**2
    assert dominant_term(4*x**5 + 6*x**3 + 2*x**2 + 5, x) == x**5
    assert dominant_term(z**3 + 3*z**2 +4*z + 8, z) == z**3
    assert dominant_term(5*z**4 + 3*z**2 + 9*z + 5, z) == z**4
    assert dominant_term(3*w**3 + 4*w**2 + 2*w + 8, w) == w**3
