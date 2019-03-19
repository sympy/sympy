from sympy.algebras.quaternion import Quaternion
from sympy import Symbol, symbols, re, im, Add, Mul, I, Abs
from sympy import cos, sin, sqrt, conjugate, log, acos, E, pi
from sympy.utilities.pytest import raises
from sympy import Matrix
from sympy import diff, integrate, trigsimp
from sympy import S, Rational
from sympy.vector.vector import Vector

x, y, z, w = symbols("x y z w")

def test_quaternion_construction():
    q = Quaternion(x, y, z, w)
    assert q + q == Quaternion(2*x, 2*y, 2*z, 2*w)

    q2 = Quaternion.from_axis_angle((sqrt(3)/3, sqrt(3)/3, sqrt(3)/3), 2*pi/3)
    assert q2 == Quaternion(Rational(1, 2), Rational(1, 2),
                            Rational(1, 2), Rational(1, 2))

    M = Matrix([[cos(x), -sin(x), 0], [sin(x), cos(x), 0], [0, 0, 1]])
    q3 = trigsimp(Quaternion.from_rotation_matrix(M))
    assert q3 == Quaternion(sqrt(2)*sqrt(cos(x) + 1)/2, 0, 0, sqrt(-2*cos(x) + 2)/2)


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

    q1 = Quaternion(3 + 4*I, 2 + 5*I, 0, 7 + 8*I, real_field = False)
    q2 = Quaternion(1, 4, 7, 8)

    assert q1 + (2 + 3*I) == Quaternion(5 + 7*I, 2 + 5*I, 0, 7 + 8*I)
    assert q2 + (2 + 3*I) == Quaternion(3, 7, 7, 8)
    assert q1 * (2 + 3*I) == \
    Quaternion((2 + 3*I)*(3 + 4*I), (2 + 3*I)*(2 + 5*I), 0, (2 + 3*I)*(7 + 8*I))
    assert q2 * (2 + 3*I) == Quaternion(-10, 11, 38, -5)


def test_quaternion_functions():
    q = Quaternion(x, y, z, w)
    q1 = Quaternion(1, 2, 3, 4)
    q0 = Quaternion(0, 0, 0, 0)

    assert conjugate(q) == Quaternion(x, -y, -z, -w)
    assert q.norm() == sqrt(w**2 + x**2 + y**2 + z**2)
    assert q.normalize() == Quaternion(x, y, z, w) / sqrt(w**2 + x**2 + y**2 + z**2)
    assert q.inverse() == Quaternion(x, -y, -z, -w) / (w**2 + x**2 + y**2 + z**2)
    raises(ValueError, lambda: q0.inverse())
    assert q.pow(2) == Quaternion(-w**2 + x**2 - y**2 - z**2, 2*x*y, 2*x*z, 2*w*x)
    assert q1.pow(-2) == Quaternion(-S(7)/225, -S(1)/225, -S(1)/150, -S(2)/225)
    assert q1.pow(-0.5) == NotImplemented

    assert q1.exp() == \
    Quaternion(E * cos(sqrt(29)),
               2 * sqrt(29) * E * sin(sqrt(29)) / 29,
               3 * sqrt(29) * E * sin(sqrt(29)) / 29,
               4 * sqrt(29) * E * sin(sqrt(29)) / 29)
    assert q1._ln() == \
    Quaternion(log(sqrt(30)),
               2 * sqrt(29) * acos(sqrt(30)/30) / 29,
               3 * sqrt(29) * acos(sqrt(30)/30) / 29,
               4 * sqrt(29) * acos(sqrt(30)/30) / 29)

    assert q1.pow_cos_sin(2) == \
    Quaternion(30 * cos(2 * acos(sqrt(30)/30)),
               60 * sqrt(29) * sin(2 * acos(sqrt(30)/30)) / 29,
               90 * sqrt(29) * sin(2 * acos(sqrt(30)/30)) / 29,
               120 * sqrt(29) * sin(2 * acos(sqrt(30)/30)) / 29)

    assert diff(Quaternion(x, x, x, x), x) == Quaternion(1, 1, 1, 1)

    assert integrate(Quaternion(x, x, x, x), x) == \
    Quaternion(x**2 / 2, x**2 / 2, x**2 / 2, x**2 / 2)

    assert Quaternion.rotate_point((1, 1, 1), q1) == (S(1) / 5, 1, S(7) / 5)


def test_quaternion_conversions():
    q1 = Quaternion(1, 2, 3, 4)

    assert q1.to_axis_angle() == ((2 * sqrt(29)/29,
                                   3 * sqrt(29)/29,
                                   4 * sqrt(29)/29),
                                   2 * acos(sqrt(30)/30))

    assert q1.to_rotation_matrix() == Matrix([[-S(2)/3, S(2)/15, S(11)/15],
                                     [S(2)/3, -S(1)/3, S(2)/3],
                                     [S(1)/3, S(14)/15, S(2)/15]])

    assert q1.to_rotation_matrix((1, 1, 1)) == Matrix([[-S(2)/3, S(2)/15, S(11)/15, S(4)/5],
                                                  [S(2)/3, -S(1)/3, S(2)/3, S(0)],
                                                       [S(1)/3, S(14)/15, S(2)/15, -S(2)/5],
                                                  [S(0), S(0), S(0), S(1)]])

    theta = symbols("theta", real=True)
    q2 = Quaternion(cos(theta/2), 0, 0, sin(theta/2))

    assert trigsimp(q2.to_rotation_matrix()) == Matrix([
                                               [cos(theta), -sin(theta), 0],
                                               [sin(theta),  cos(theta), 0],
                                               [0,           0,          1]])

    assert q2.to_axis_angle() == ((0, 0, sin(theta/2)/Abs(sin(theta/2))),
                                   2*acos(cos(theta/2)))

    assert trigsimp(q2.to_rotation_matrix((1, 1, 1))) == Matrix([
               [cos(theta), -sin(theta), 0, sin(theta) - cos(theta) + 1],
               [sin(theta),  cos(theta), 0, -sin(theta) - cos(theta) + 1],
               [0,           0,          1,  0],
               [0,           0,          0,  1]])


def test_quaternion_rotation_iss1593():
    """
    There was a sign mistake in the definition,
    of the rotation matrix. This tests that particular sign mistake.
    See issue 1593 for reference.
    See wikipedia
    https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
    for the correct definition
    """
    q = Quaternion(cos(x/2), sin(x/2), 0, 0)
    assert(trigsimp(q.to_rotation_matrix()) == Matrix([
                [1,      0,      0],
                [0, cos(x), -sin(x)],
                [0, sin(x), cos(x)]]))


def test_quaternion_differentiation():
    q = Quaternion(1, 4, 5, 7)
    v1 = Vector(1, 2, 3)
    v2 = Vector(1, 4, 7)
    v3 = Vector(9, 5, 8)
    v4 = Vector(3, 5, 2)
    assert q.differentiate(v1) == Quaternion(-17.5, 1.0, -1.5, 3.0)
    assert q.differentiate(v2) == Quaternion(-36.5, 4.0, -8.5, 9.0)
    assert q.differentiate(v3) == Quaternion(-58.5, 7.0, 18.0, -8.5)
    assert q.differentiate(v4) == Quaternion(-25.5, -11.0, 9.0, 3.5)


def test_quaternion_integration():
    x = Symbol('x')
    q1 = Quaternion(2, 9, 6, 3)
    q2 = Quaternion(x-1, x-2, x-3, x-4)
    q3 = Quaternion(1/x, 1/x**2, 1/x**3, 1/x**4)
    q4 = Quaternion(sin(x), sin(2*x), cos(x), cos(2*x))

    assert q1.integrate(x) == Quaternion(2*x, 9*x, 6*x, 3*x)
    assert q2.integrate(x) == Quaternion(x**2/2-x, x**2/2 -2*x, x**2/2 -3*x, x**2/2 - 4*x)
    assert q3.integrate(x) == Quaternion(log(x), -1/x, -1/(2*x**2), -1/(3*x**3))
    assert q4.integrate(x) == Quaternion(-cos(x), -cos(2*x)/2, sin(x), sin(2*x)/2)
