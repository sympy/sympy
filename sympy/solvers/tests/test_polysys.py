"""Tests for solvers of systems of polynomial equations. """

from sympy import S, symbols, Integer, Rational, sqrt, I, Poly, QQ, \
                  all, flatten

from sympy.abc import x, y, z
from sympy.solvers.polysys import solve_poly_system, solve_triangulated
from sympy.utilities.pytest import raises


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

    f_1 = x**2 + y + z - 1
    f_2 = x + y**2 + z - 1
    f_3 = x + y + z**2 - 1

    a, b = -sqrt(2) - 1, sqrt(2) - 1

    assert solve_poly_system([f_1, f_2, f_3], x, y, z) == \
        [(a, a, a), (0, 0, 1), (0, 1, 0), (b, b, b), (1, 0, 0)]

    solution = [(1, -1), (1, 1)]

    assert solve_poly_system([Poly(x**2 - y**2), Poly(x - 1)]) == solution
    assert solve_poly_system([x**2 - y**2, x - 1], x, y) == solution
    assert solve_poly_system([x**2 - y**2, x - 1]) == solution

    raises(NotImplementedError, "solve_poly_system([x**3-y**3], x, y)")

def test_solve_biquadratic():
    x0, y0, x1, y1, r = symbols('x0 y0 x1 y1 r')

    f_1 = (x - 1)**2 + (y - 1)**2 - r**2
    f_2 = (x - 2)**2 + (y - 2)**2 - r**2

    assert solve_poly_system([f_1, f_2], x, y) == \
        [(S(3)/2 + (-1 + 2*r**2)**(S(1)/2)/2, S(3)/2 - (-1 + 2*r**2)**(S(1)/2)/2),
         (S(3)/2 - (-1 + 2*r**2)**(S(1)/2)/2, S(3)/2 + (-1 + 2*r**2)**(S(1)/2)/2)]

    f_1 = (x - 1)**2 + (y - 2)**2 - r**2
    f_2 = (x - 1)**2 + (y - 1)**2 - r**2

    assert solve_poly_system([f_1, f_2], x, y) == \
        [(1 + (-((1 - 2*r)*(1 + 2*r)))**(S(1)/2)/2, S(3)/2),
         (1 - (-((1 - 2*r)*(1 + 2*r)))**(S(1)/2)/2, S(3)/2)]

    query = lambda expr: expr.is_Pow and expr.exp is S.Half

    f_1 = (x - 1 )**2 + (y - 2)**2 - r**2
    f_2 = (x - x1)**2 + (y - 1)**2 - r**2

    result = solve_poly_system([f_1, f_2], x, y)

    assert len(result) == 2 and all(len(r) == 2 for r in result)
    assert all(len(r.find(query)) == 1 for r in flatten(result))

    f_1 = (x - x0)**2 + (y - y0)**2 - r**2
    f_2 = (x - x1)**2 + (y - y1)**2 - r**2

    result = solve_poly_system([f_1, f_2], x, y)

    assert len(result) == 2 and all(len(r) == 2 for r in result)
    assert all(len(r.find(query)) == 1 for r in flatten(result))

def test_solve_triangualted():
    f_1 = x**2 + y + z - 1
    f_2 = x + y**2 + z - 1
    f_3 = x + y + z**2 - 1

    a, b = -sqrt(2) - 1, sqrt(2) - 1

    assert solve_triangulated([f_1, f_2, f_3], x, y, z) == \
        [(0, 0, 1), (0, 1, 0), (1, 0, 0)]

    dom = QQ.algebraic_field(sqrt(2))

    assert solve_triangulated([f_1, f_2, f_3], x, y, z, domain=dom) == \
        [(a, a, a), (0, 0, 1), (0, 1, 0), (b, b, b), (1, 0, 0)]
