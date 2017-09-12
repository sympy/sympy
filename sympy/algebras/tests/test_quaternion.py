from sympy.algebras.quaternion import Quaternion
from sympy import symbols, re, im, Add, Mul, I
from sympy import cos, sin, sqrt, conjugate, exp, log, acos, E
from sympy import Matrix
from sympy import diff, integrate


x, y, z, w = symbols("x y z w")


def test_quaternion_construction():
    q = Quaternion(x, y, z, w)
    assert q + q == Quaternion(2 * x, 2 * y, 2 * z, 2 * w)

    q2 = Quaternion.from_axis_angle((1, 4, 5), 2)
    assert q2 == Quaternion(
        cos(1),
        sqrt(42) * sin(1) / 42,
        2 * sqrt(42) * sin(1) / 21,
        5 * sqrt(42) * sin(1) / 42)

    q3 = Quaternion.from_matrix(Matrix([[1, 2, 3], [0, 1, 4], [5, 6, 0]]))
    assert q3 == Quaternion(sqrt(3) / 2, 1 / 2, -1 / 2, 0)

    q4 = Quaternion(3 + 4 * I, 0, 0, 0)
    assert q4 is None


def test_quaternion_complex_real_addition():
    a = symbols("a", complex=True)
    b = symbols("b", real=True)
    # This symbol is not complex:
    c = symbols("c", commutative=False)

    q = Quaternion(x, y, z, w)
    assert a + q == Quaternion(x + re(a), y + im(a), z, w)
    assert 1 + q == Quaternion(1 + x, y, z, w)
    assert I + q == Quaternion(x, 1 + y, z, w)
    assert b + q == Quaternion(x + b, y, z, w)
    assert c + q == Add(c, Quaternion(x, y, z, w), evaluate=False)
    assert q * c == Mul(Quaternion(x, y, z, w), c, evaluate=False)
    assert c * q == Mul(c, Quaternion(x, y, z, w), evaluate=False)

    assert -q == Quaternion(-x, -y, -z, -w)


def test_quaternion_functions():
    q = Quaternion(x, y, z, w)
    q1 = Quaternion(1, 2, 3, 4)

    assert conjugate(q) == Quaternion(x, -y, -z, -w)
    assert q.norm() == sqrt(w**2 + x**2 + y**2 + z**2)
    assert q.normalize() == Quaternion(x, y, z, w) / \
        sqrt(w**2 + x**2 + y**2 + z**2)
    assert q.inverse() == Quaternion(x, -y, -z, -w) / (w**2 + x**2 + y**2 + z**2)

    assert q1.exp() == Quaternion(E * cos(sqrt(29)),
                                  2 * sqrt(29) * E * sin(sqrt(29)) / 29,
                                  3 * sqrt(29) * E * sin(sqrt(29)) / 29,
                                  4 * sqrt(29) * E * sin(sqrt(29)) / 29)

    assert q1.ln() == Quaternion(log(sqrt(30)),
                                 2 * sqrt(29) * acos(sqrt(30) / 30) / 29,
                                 3 * sqrt(29) * acos(sqrt(30) / 30) / 29,
                                 4 * sqrt(29) * acos(sqrt(30) / 30) / 29)

    assert q1.pow_cos_sin(2) == Quaternion(30 * cos(2 * acos(sqrt(30) / 30)),
                                           60 * sqrt(29) * sin(2 * acos(sqrt(30) / 30)) / 29,
                                           90 * sqrt(29) * sin(2 * acos(sqrt(30) / 30)) / 29,
                                           120 * sqrt(29) * sin(2 * acos(sqrt(30) / 30)) / 29)

    assert diff(Quaternion(x, x, x, x), x) == Quaternion(1, 1, 1, 1)
    assert integrate(Quaternion(x,x,x,x),x) == Quaternion(x**2 / 2, x**2 / 2, x**2 / 2, x**2 / 2)

    assert Quaternion.rotate((1, 1, 1), q1) == (1 / 5, 1, 7 / 5)


def test_quaternion_conversions():
    q1 = Quaternion(1, 2, 3, 4)

    assert q1.to_axis_angle() == ((2 * sqrt(29) / 29, 3 * sqrt(29) / 29,
                                   4 * sqrt(29) / 29), 2 * acos(sqrt(30) / 30))

    assert q1.to_matrix() == Matrix([[-2 / 5, 2 / 5, 11 / 5],
                                     [2, -1 / 5, 2],
                                     [1, 14 / 5, 0]])

    assert q1.to_matrix_4x4((1, 1, 1)) == Matrix([[-2 / 5, 2 / 5, 11 / 5, -6 / 5],
                                                  [2, -1 / 5, 2, -14 / 5],
                                                  [1, 14 / 5, 0, -14 / 5],
                                                  [0, 0, 0, 1]])
