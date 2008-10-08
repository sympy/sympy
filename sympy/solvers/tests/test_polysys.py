
from sympy import S, symbols, Integer, Rational, Poly, \
    PolynomialError, sqrt, I, solve_poly_system, raises

x, y, z = symbols('xyz')

def test_solve_poly_system():
    assert solve_poly_system([x-1], x) == [(S.One,)]

    assert solve_poly_system([y - x, y - x - 1], x, y) == None

    assert solve_poly_system([y - x**2, y + x**2], x, y) == [(S.Zero, S.Zero)]

    assert solve_poly_system([2*x - 3, 3*y/2 - 2*x, z - 5*y], x, y, z) == \
        [(Rational(3, 2), Integer(2), Integer(10))]

    assert solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y) == \
       [(0, 0), (2, -2**S.Half), (2, 2**S.Half)]

    assert solve_poly_system([y - x**2, y + x**2 + 1], x, y) == \
           [(I*S.Half**S.Half, -S.Half), (-I*S.Half**S.Half, -S.Half)]

    a, b = -1 - 2**S.Half, -1 + 2**S.Half

    assert solve_poly_system([x**2+y+z-1, x+y**2+z-1, x+y+z**2-1], x, y, z) == \
        [(a, a, a), (0, 0, 1), (0, 1, 0), (b, b, b), (1, 0, 0)]

    raises(PolynomialError, "solve_poly_system([x**3-y**3], x, y)")
