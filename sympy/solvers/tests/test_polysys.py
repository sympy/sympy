from sympy import S, symbols, Integer, Rational, sqrt, I, solve_poly_system, Poly

from sympy.utilities.pytest import raises

from sympy.abc import x, y, z

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

    solution = [(1, -1), (1, 1)]

    assert solve_poly_system([Poly(x**2 - y**2), Poly(x - 1)]) == solution
    assert solve_poly_system([x**2 - y**2, x - 1], x, y) == solution
    assert solve_poly_system([x**2 - y**2, x - 1]) == solution

    raises(NotImplementedError, "solve_poly_system([x**3-y**3], x, y)")

def test_solve_poly_system_slow():
    f_1 = x**2 + y + z - 1
    f_2 = x + y**2 + z - 1
    f_3 = x + y + z**2 - 1

    a, b = -sqrt(2) - 1, sqrt(2) - 1

    assert solve_poly_system([f_1, f_2, f_3], x, y, z) == \
        [(a, a, a), (0, 0, 1), (0, 1, 0), (b, b, b), (1, 0, 0)]
