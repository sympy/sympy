from sympy import sin, cos
from sympy import pi
from sympy.series.kauers import Delta
from sympy.abc import x
from sympy import expand

def test_Delta():
    assert Delta(x**2 + 2*x + 1) == 2*x + 3
    assert Delta(x**3 + 2*x**2 + 3*x +5) == 3*x**2 + 7*x + 6
    assert Delta(x**2 - 2*x + 3) == 2*x - 1
    assert Delta(x**2 + 3*x -2) == 2*x + 4
    assert Delta(sin(x), pi/6) == -sin(x) + sin(x + pi/6) 
    assert Delta(cos(x), pi/3)
    assert Delta(x**2 - 2*x + 3, 2) == 4*x
    assert Delta(x**2 - 2*x + 3, 3) == 6*x + 3
    

