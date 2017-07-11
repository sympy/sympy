from sympy import Curve, E, line_integrate, ln, sqrt, symbols

s, t, x, y, z = symbols('s,t,x,y,z')


def test_lineintegral():
    c = Curve([E**t + 1, E**t - 1], (t, 0, ln(2)))
    assert line_integrate(x + y, c, [x, y]) == 3*sqrt(2)
