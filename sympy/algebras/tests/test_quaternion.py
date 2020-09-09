from sympy import symbols, re, im, sign, I, Abs, Symbol, \
     cos, sin, sqrt, conjugate, log, acos, E, pi, \
     Matrix, diff, integrate, trigsimp, S, Rational
from sympy.algebras.quaternion import Quaternion
from sympy.testing.pytest import raises

w, x, y, z = symbols('w:z')
phi = symbols('phi')

def test_quaternion_construction():
    q = Quaternion(w, x, y, z)
    assert q + q == Quaternion(2*w, 2*x, 2*y, 2*z)

    q2 = Quaternion.from_axis_angle((sqrt(3)/3, sqrt(3)/3, sqrt(3)/3),
                                    pi*Rational(2, 3))
    assert q2 == Quaternion(S.Half, S.Half,
                            S.Half, S.Half)

    M = Matrix([[cos(phi), -sin(phi), 0], [sin(phi), cos(phi), 0], [0, 0, 1]])
    q3 = trigsimp(Quaternion.from_rotation_matrix(M))
    assert q3 == Quaternion(sqrt(2)*sqrt(cos(phi) + 1)/2, 0, 0, sqrt(2 - 2*cos(phi))*sign(sin(phi))/2)

    nc = Symbol('nc', commutative=False)
    raises(ValueError, lambda: Quaternion(w, x, nc, z))


def test_quaternion_complex_real_addition():
    a = symbols("a", complex=True)
    b = symbols("b", real=True)
    # This symbol is not complex:
    c = symbols("c", commutative=False)

    q = Quaternion(w, x, y, z)
    assert a + q == Quaternion(w + re(a), x + im(a), y, z)
    assert 1 + q == Quaternion(1 + w, x, y, z)
    assert I + q == Quaternion(w, 1 + x, y, z)
    assert b + q == Quaternion(w + b, x, y, z)
    raises(ValueError, lambda: c + q)
    raises(ValueError, lambda: q * c)
    raises(ValueError, lambda: c * q)

    assert -q == Quaternion(-w, -x, -y, -z)

    q1 = Quaternion(3 + 4*I, 2 + 5*I, 0, 7 + 8*I, real_field = False)
    q2 = Quaternion(1, 4, 7, 8)

    assert q1 + (2 + 3*I) == Quaternion(5 + 7*I, 2 + 5*I, 0, 7 + 8*I)
    assert q2 + (2 + 3*I) == Quaternion(3, 7, 7, 8)
    assert q1 * (2 + 3*I) == \
    Quaternion((2 + 3*I)*(3 + 4*I), (2 + 3*I)*(2 + 5*I), 0, (2 + 3*I)*(7 + 8*I))
    assert q2 * (2 + 3*I) == Quaternion(-10, 11, 38, -5)

    q1 = Quaternion(1, 2, 3, 4)
    q0 = Quaternion(0, 0, 0, 0)
    assert q1 + q0 == q1
    assert q1 - q0 == q1
    assert q1 - q1 == q0


def test_quaternion_functions():
    q = Quaternion(w, x, y, z)
    q1 = Quaternion(1, 2, 3, 4)
    q0 = Quaternion(0, 0, 0, 0)

    assert conjugate(q) == Quaternion(w, -x, -y, -z)
    assert q.norm() == sqrt(w**2 + x**2 + y**2 + z**2)
    assert q.normalize() == Quaternion(w, x, y, z) / sqrt(w**2 + x**2 + y**2 + z**2)
    assert q.inverse() == Quaternion(w, -x, -y, -z) / (w**2 + x**2 + y**2 + z**2)
    assert q.inverse() == q.pow(-1)
    raises(ValueError, lambda: q0.inverse())
    assert q.pow(2) == Quaternion(w**2 - x**2 - y**2 - z**2, 2*w*x, 2*w*y, 2*w*z)
    assert q**(2) == Quaternion(w**2 - x**2 - y**2 - z**2, 2*w*x, 2*w*y, 2*w*z)
    assert q1.pow(-2) == Quaternion(Rational(-7, 225), Rational(-1, 225), Rational(-1, 150), Rational(-2, 225))
    assert q1**(-2) == Quaternion(Rational(-7, 225), Rational(-1, 225), Rational(-1, 150), Rational(-2, 225))
    assert q1.pow(-0.5) == NotImplemented
    raises(TypeError, lambda: q1**(-0.5))

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

    assert Quaternion.rotate_point((1, 1, 1), q1) == (S.One / 5, 1, S(7) / 5)
    n = Symbol('n')
    raises(TypeError, lambda: q1**n)
    n = Symbol('n', integer=True)
    raises(TypeError, lambda: q1**n)


def test_quaternion_conversions():
    q1 = Quaternion(1, 2, 3, 4)

    assert q1.to_axis_angle() == ((2 * sqrt(29)/29,
                                   3 * sqrt(29)/29,
                                   4 * sqrt(29)/29),
                                   2 * acos(sqrt(30)/30))

    assert q1.to_rotation_matrix() == Matrix([[Rational(-2, 3), Rational(2, 15), Rational(11, 15)],
                                              [Rational(2, 3), Rational(-1, 3), Rational(2, 3)],
                                              [Rational(1, 3), Rational(14, 15), Rational(2, 15)]])

    assert q1.to_rotation_matrix((1, 1, 1)) == Matrix([[Rational(-2, 3), Rational(2, 15), Rational(11, 15), Rational(4, 5)],
                                                       [Rational(2, 3), Rational(-1, 3), Rational(2, 3), S.Zero],
                                                       [Rational(1, 3), Rational(14, 15), Rational(2, 15), Rational(-2, 5)],
                                                       [S.Zero, S.Zero, S.Zero, S.One]])

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
    q = Quaternion(cos(phi/2), sin(phi/2), 0, 0)
    assert(trigsimp(q.to_rotation_matrix()) == Matrix([
                [1,        0,         0],
                [0, cos(phi), -sin(phi)],
                [0, sin(phi),  cos(phi)]]))


def test_quaternion_multiplication():
    q1 = Quaternion(3 + 4*I, 2 + 5*I, 0, 7 + 8*I, real_field = False)
    q2 = Quaternion(1, 2, 3, 5)
    q3 = Quaternion(1, 1, 1, y)

    assert Quaternion._generic_mul(4, 1) == 4
    assert Quaternion._generic_mul(4, q1) == Quaternion(12 + 16*I, 8 + 20*I, 0, 28 + 32*I)
    assert q2.mul(2) == Quaternion(2, 4, 6, 10)
    assert q2.mul(q3) == Quaternion(-5*y - 4, 3*y - 2, 9 - 2*y, y + 4)
    assert q2.mul(q3) == q2*q3

    z = symbols('z', complex=True)
    z_quat = Quaternion(re(z), im(z), 0, 0)
    q = Quaternion(*symbols('q:4', real=True))

    assert z * q == z_quat * q
    assert q * z == q * z_quat


def test_issue_16318():
    #for rtruediv
    q0 = Quaternion(0, 0, 0, 0)
    raises(ValueError, lambda: 1/q0)
    #for rotate_point
    q = Quaternion(1, 2, 3, 4)
    (axis, angle) = q.to_axis_angle()
    assert Quaternion.rotate_point((1, 1, 1), (axis, angle)) == (S.One / 5, 1, S(7) / 5)
    #test for to_axis_angle
    q = Quaternion(-1, 1, 1, 1)
    axis = (-sqrt(3)/3, -sqrt(3)/3, -sqrt(3)/3)
    angle = 2*pi/3
    assert (axis, angle) == q.to_axis_angle()
