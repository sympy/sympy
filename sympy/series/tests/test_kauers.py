from sympy.series.kauers import Delta
from sympy.abc import x, y, z, w, n
from sympy import  sin, cos
from sympy import  pi


def test_Delta():
    assert Delta(x**2 + 2*x + 1, x) == 2*x + 3
    assert Delta(y**3 + 2*y**2 + 3*y +5, y) == 3*y**2 + 7*y + 6
    assert Delta(z**2 - 2*z + 3, z) == 2*z - 1
    assert Delta(w**2 + 3*w -2, w) == 2*w + 4
    assert Delta(sin(x), x,  pi/6) == -sin(x) + sin(x + pi/6) 
    assert Delta(cos(y),y,  pi/3)
    assert Delta(x**2 - 2*x + 3, x,  2) == 4*x
    assert Delta(n**2 - 2*n + 3, n, 3) == 6*n + 3
    
