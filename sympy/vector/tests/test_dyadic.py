from sympy import sin, cos, symbols, pi, ImmutableMatrix as Matrix
from sympy.vector import (CoordSysCartesian, Vector, Dyadic)


A = CoordSysCartesian('A')


def test_dyadic():
    d1 = A.i | A.i
    d2 = A.j | A.j
    d3 = A.i | A.j
    assert d1 * 0 == 0
    assert d1 != 0
    assert d1 * 2 == 2 * (A.i | A.i)
    assert d1 / 2. == 0.5 * d1
    assert d1 & (0 * d1) == 0
    assert d1 & d2 == 0
    assert d1 & A.i == A.i
    assert d1 ^ A.i == 0
    assert d1 ^ A.j == A.i | A.k
    assert d1 ^ A.k == - A.i | A.j
    assert d2 ^ A.i == - A.j | A.k
    assert A.i ^ d1 == 0
    assert A.j ^ d1 == - A.k | A.i
    assert A.k ^ d1 == A.j | A.i
    assert A.i & d1 == A.i
    assert A.j & d1 == 0
    assert A.j & d2 == A.j
    assert d1 & d3 == A.i | A.j
    assert d3 & d1 == 0
    assert d1.dt(A) == 0
    q = symbols('q')
    B = A.orient_new('B', 'Axis', [q, A.k])
    assert d1.express(B) == d1.express(B, B)
    assert d1.express(B) == ((cos(q)**2) * (B.i | B.i) + (-sin(q) * cos(q)) *
            (B.i | B.j) + (-sin(q) * cos(q)) * (B.j | B.i) + (sin(q)**2) *
            (B.j | B.j))
    assert express(d1, B, A) == (cos(q)) * (B.i | A.i) + (-sin(q)) * (B.j | A.i)
    assert express(d1, A, B) == (cos(q)) * (A.i | B.i) + (-sin(q)) * (A.i | B.j)
    assert d1.to_matrix(A) == Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert d1.to_matrix(A, B) == Matrix([[cos(q), -sin(q), 0],
                                         [0, 0, 0],
                                         [0, 0, 0]])
    assert d3.to_matrix(A) == Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    a, b, c, d, e, f = symbols('a, b, c, d, e, f')
    v1 = a * A.i + b * A.j + c * A.k
    v2 = d * A.i + e * A.j + f * A.k
    d4 = v1.outer(v2)
    assert d4.to_matrix(A) == Matrix([[a * d, a * e, a * f],
                                      [b * d, b * e, b * f],
                                      [c * d, c * e, c * f]])
    d5 = v1.outer(v1)
    C = A.orient_new('C', 'Axis', [q, A.i])
    for expected, actual in zip(C.rotation_matrix(A) * d5.to_matrix(A) * \
                               C.rotation_matrix(A).T, d5.to_matrix(C)):
        assert (expected - actual).simplify() == 0


def test_dyadic_simplify():
    x, y, z, k, n, m, w, f, s, A = symbols('x, y, z, k, n, m, w, f, s, A')
    N = CoordSysCartesian('N')

    dy = N.i | N.i
    test1 = (1 / x + 1 / y) * dy
    assert (N.i & test1 & N.i) != (x + y) / (x * y)
    test1 = test1.simplify()
    assert (N.i & test1 & N.i) == (x + y) / (x * y)

    test2 = (A**2 * s**4 / (4 * pi * k * m**3)) * dy
    test2 = test2.simplify()
    assert (N.i & test2 & N.i) == (A**2 * s**4 / (4 * pi * k * m**3))

    test3 = ((4 + 4 * x - 2 * (2 + 2 * x)) / (2 + 2 * x)) * dy
    test3 = test3.simplify()
    assert (N.i & test3 & N.i) == 0

    test4 = ((-4 * x * y**2 - 2 * y**3 - 2 * x**2 * y) / (x + y)**2) * dy
    test4 = test4.simplify()
    assert (N.i & test4 & N.i) == -2 * y
