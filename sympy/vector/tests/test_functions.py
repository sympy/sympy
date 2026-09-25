from __future__ import annotations
from sympy.vector.vector import Vector
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.functions import (
    express, matrix_to_vector, matrix_to_dyadic, orthogonalize)
from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.matrices.exceptions import ShapeError
from sympy.matrices.immutable import ImmutableDenseMatrix as Matrix
from sympy.testing.pytest import raises

N = CoordSys3D('N')
q1, q2, q3, q4, q5 = symbols('q1 q2 q3 q4 q5')
A = N.orient_new_axis('A', q1, N.k)  # type: ignore
B = A.orient_new_axis('B', q2, A.i)
C = B.orient_new_axis('C', q3, B.j)


def test_express():
    assert express(Vector.zero, N) == Vector.zero
    assert express(S.Zero, N) is S.Zero
    assert express(A.i, C) == cos(q3)*C.i + sin(q3)*C.k
    assert express(A.j, C) == sin(q2)*sin(q3)*C.i + cos(q2)*C.j - \
        sin(q2)*cos(q3)*C.k
    assert express(A.k, C) == -sin(q3)*cos(q2)*C.i + sin(q2)*C.j + \
        cos(q2)*cos(q3)*C.k
    assert express(A.i, N) == cos(q1)*N.i + sin(q1)*N.j
    assert express(A.j, N) == -sin(q1)*N.i + cos(q1)*N.j
    assert express(A.k, N) == N.k
    assert express(A.i, A) == A.i
    assert express(A.j, A) == A.j
    assert express(A.k, A) == A.k
    assert express(A.i, B) == B.i
    assert express(A.j, B) == cos(q2)*B.j - sin(q2)*B.k
    assert express(A.k, B) == sin(q2)*B.j + cos(q2)*B.k
    assert express(A.i, C) == cos(q3)*C.i + sin(q3)*C.k
    assert express(A.j, C) == sin(q2)*sin(q3)*C.i + cos(q2)*C.j - \
        sin(q2)*cos(q3)*C.k
    assert express(A.k, C) == -sin(q3)*cos(q2)*C.i + sin(q2)*C.j + \
        cos(q2)*cos(q3)*C.k
    # Check to make sure UnitVectors get converted properly
    assert express(N.i, N) == N.i
    assert express(N.j, N) == N.j
    assert express(N.k, N) == N.k
    assert express(N.i, A) == (cos(q1)*A.i - sin(q1)*A.j)
    assert express(N.j, A) == (sin(q1)*A.i + cos(q1)*A.j)
    assert express(N.k, A) == A.k
    assert express(N.i, B) == (cos(q1)*B.i - sin(q1)*cos(q2)*B.j +
            sin(q1)*sin(q2)*B.k)
    assert express(N.j, B) == (sin(q1)*B.i + cos(q1)*cos(q2)*B.j -
            sin(q2)*cos(q1)*B.k)
    assert express(N.k, B) == (sin(q2)*B.j + cos(q2)*B.k)
    assert express(N.i, C) == (
        (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*C.i -
        sin(q1)*cos(q2)*C.j +
        (sin(q3)*cos(q1) + sin(q1)*sin(q2)*cos(q3))*C.k)
    assert express(N.j, C) == (
        (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.i +
        cos(q1)*cos(q2)*C.j +
        (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.k)
    assert express(N.k, C) == (-sin(q3)*cos(q2)*C.i + sin(q2)*C.j +
            cos(q2)*cos(q3)*C.k)

    assert express(A.i, N) == (cos(q1)*N.i + sin(q1)*N.j)
    assert express(A.j, N) == (-sin(q1)*N.i + cos(q1)*N.j)
    assert express(A.k, N) == N.k
    assert express(A.i, A) == A.i
    assert express(A.j, A) == A.j
    assert express(A.k, A) == A.k
    assert express(A.i, B) == B.i
    assert express(A.j, B) == (cos(q2)*B.j - sin(q2)*B.k)
    assert express(A.k, B) == (sin(q2)*B.j + cos(q2)*B.k)
    assert express(A.i, C) == (cos(q3)*C.i + sin(q3)*C.k)
    assert express(A.j, C) == (sin(q2)*sin(q3)*C.i + cos(q2)*C.j -
            sin(q2)*cos(q3)*C.k)
    assert express(A.k, C) == (-sin(q3)*cos(q2)*C.i + sin(q2)*C.j +
            cos(q2)*cos(q3)*C.k)

    assert express(B.i, N) == (cos(q1)*N.i + sin(q1)*N.j)
    assert express(B.j, N) == (-sin(q1)*cos(q2)*N.i +
            cos(q1)*cos(q2)*N.j + sin(q2)*N.k)
    assert express(B.k, N) == (sin(q1)*sin(q2)*N.i -
            sin(q2)*cos(q1)*N.j + cos(q2)*N.k)
    assert express(B.i, A) == A.i
    assert express(B.j, A) == (cos(q2)*A.j + sin(q2)*A.k)
    assert express(B.k, A) == (-sin(q2)*A.j + cos(q2)*A.k)
    assert express(B.i, B) == B.i
    assert express(B.j, B) == B.j
    assert express(B.k, B) == B.k
    assert express(B.i, C) == (cos(q3)*C.i + sin(q3)*C.k)
    assert express(B.j, C) == C.j
    assert express(B.k, C) == (-sin(q3)*C.i + cos(q3)*C.k)

    assert express(C.i, N) == (
        (cos(q1)*cos(q3) - sin(q1)*sin(q2)*sin(q3))*N.i +
        (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*N.j -
        sin(q3)*cos(q2)*N.k)
    assert express(C.j, N) == (
        -sin(q1)*cos(q2)*N.i + cos(q1)*cos(q2)*N.j + sin(q2)*N.k)
    assert express(C.k, N) == (
        (sin(q3)*cos(q1) + sin(q1)*sin(q2)*cos(q3))*N.i +
        (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*N.j +
        cos(q2)*cos(q3)*N.k)
    assert express(C.i, A) == (cos(q3)*A.i + sin(q2)*sin(q3)*A.j -
            sin(q3)*cos(q2)*A.k)
    assert express(C.j, A) == (cos(q2)*A.j + sin(q2)*A.k)
    assert express(C.k, A) == (sin(q3)*A.i - sin(q2)*cos(q3)*A.j +
            cos(q2)*cos(q3)*A.k)
    assert express(C.i, B) == (cos(q3)*B.i - sin(q3)*B.k)
    assert express(C.j, B) == B.j
    assert express(C.k, B) == (sin(q3)*B.i + cos(q3)*B.k)
    assert express(C.i, C) == C.i
    assert express(C.j, C) == C.j
    assert express(C.k, C) == C.k == (C.k)

    #  Check to make sure Vectors get converted back to UnitVectors
    assert N.i == express((cos(q1)*A.i - sin(q1)*A.j), N).simplify()
    assert N.j == express((sin(q1)*A.i + cos(q1)*A.j), N).simplify()
    assert N.i == express((cos(q1)*B.i - sin(q1)*cos(q2)*B.j +
            sin(q1)*sin(q2)*B.k), N).simplify()
    assert N.j == express((sin(q1)*B.i + cos(q1)*cos(q2)*B.j -
        sin(q2)*cos(q1)*B.k), N).simplify()
    assert N.k == express((sin(q2)*B.j + cos(q2)*B.k), N).simplify()


    assert A.i == express((cos(q1)*N.i + sin(q1)*N.j), A).simplify()
    assert A.j == express((-sin(q1)*N.i + cos(q1)*N.j), A).simplify()

    assert A.j == express((cos(q2)*B.j - sin(q2)*B.k), A).simplify()
    assert A.k == express((sin(q2)*B.j + cos(q2)*B.k), A).simplify()

    assert A.i == express((cos(q3)*C.i + sin(q3)*C.k), A).simplify()
    assert A.j == express((sin(q2)*sin(q3)*C.i + cos(q2)*C.j -
            sin(q2)*cos(q3)*C.k), A).simplify()

    assert A.k == express((-sin(q3)*cos(q2)*C.i + sin(q2)*C.j +
            cos(q2)*cos(q3)*C.k), A).simplify()
    assert B.i == express((cos(q1)*N.i + sin(q1)*N.j), B).simplify()
    assert B.j == express((-sin(q1)*cos(q2)*N.i +
            cos(q1)*cos(q2)*N.j + sin(q2)*N.k), B).simplify()

    assert B.k == express((sin(q1)*sin(q2)*N.i -
            sin(q2)*cos(q1)*N.j + cos(q2)*N.k), B).simplify()

    assert B.j == express((cos(q2)*A.j + sin(q2)*A.k), B).simplify()
    assert B.k == express((-sin(q2)*A.j + cos(q2)*A.k), B).simplify()
    assert B.i == express((cos(q3)*C.i + sin(q3)*C.k), B).simplify()
    assert B.k == express((-sin(q3)*C.i + cos(q3)*C.k), B).simplify()
    assert C.i == express((cos(q3)*A.i + sin(q2)*sin(q3)*A.j -
            sin(q3)*cos(q2)*A.k), C).simplify()
    assert C.j == express((cos(q2)*A.j + sin(q2)*A.k), C).simplify()
    assert C.k == express((sin(q3)*A.i - sin(q2)*cos(q3)*A.j +
            cos(q2)*cos(q3)*A.k), C).simplify()
    assert C.i == express((cos(q3)*B.i - sin(q3)*B.k), C).simplify()
    assert C.k == express((sin(q3)*B.i + cos(q3)*B.k), C).simplify()


def test_matrix_to_vector():
    m = Matrix([[1], [2], [3]])
    assert matrix_to_vector(m, C) == C.i + 2*C.j + 3*C.k
    m = Matrix([[0], [0], [0]])
    assert matrix_to_vector(m, N) == matrix_to_vector(m, C) == \
           Vector.zero
    m = Matrix([[q1], [q2], [q3]])
    assert matrix_to_vector(m, N) == q1*N.i + q2*N.j + q3*N.k


def test_matrix_to_dyadic():
    raises(ShapeError, lambda: matrix_to_dyadic(Matrix([0]), N))
    raises(
        ShapeError,
        lambda: matrix_to_dyadic(Matrix([[0, 1, 2], [3, 4, 5]]), N))

    m = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    expected = (
        1 * (N.i | N.i) + 2 * (N.i | N.j) + 3 * (N.i | N.k) +
        4 * (N.j | N.i) + 5 * (N.j | N.j) + 6 * (N.j | N.k) +
        7 * (N.k | N.i) + 8 * (N.k | N.j) + 9 * (N.k | N.k)
    )
    assert matrix_to_dyadic(m, N) == expected


def test_orthogonalize():
    C = CoordSys3D('C')
    a, b = symbols('a b', integer=True)
    i, j, k = C.base_vectors()
    v1 = i + 2*j
    v2 = 2*i + 3*j
    v3 = 3*i + 5*j
    v4 = 3*i + j
    v5 = 2*i + 2*j
    v6 = a*i + b*j
    v7 = 4*a*i + 4*b*j
    assert orthogonalize(v1, v2) == [C.i + 2*C.j, C.i*Rational(2, 5) + -C.j/5]
    # from wikipedia
    assert orthogonalize(v4, v5, orthonormal=True) == \
        [(3*sqrt(10))*C.i/10 + (sqrt(10))*C.j/10, (-sqrt(10))*C.i/10 + (3*sqrt(10))*C.j/10]
    raises(ValueError, lambda: orthogonalize(v1, v2, v3))
    raises(ValueError, lambda: orthogonalize(v6, v7))


def test_express_from_spherical_to_cartesian():
    C = CoordSys3D("C")
    S = C.create_new("S", transformation="spherical")
    r, theta, phi = S.base_scalars()
    a, b, c = symbols("a:c")

    # vector along the radial direction
    v1 = 1 * S.i
    assert express(v1, C) == (
        (sin(theta) * cos(phi)) * C.i +
        (sin(theta) * sin(phi)) * C.j +
        cos(theta) * C.k)
    # vector along the polar direction
    v2 = 1 * S.j
    assert express(v2, C) == (
        (cos(theta) * cos(phi)) * C.i +
        (cos(theta) * sin(phi)) * C.j +
        (-sin(theta)) * C.k)
    # vector along the azimuthal direction
    v3 = 1 * S.k
    assert express(v3, C) == -sin(phi) * C.i + cos(phi) * C.j
    # arbitrary vector expressed in spherical coordinates
    v4 = a * S.i + b * S.j + c * S.k
    assert express(v4, C) == (
        (a*cos(phi)*sin(theta) + b*cos(phi)*cos(theta) - c*sin(phi)) * C.i +
        (a*sin(phi)*sin(theta) + b*sin(phi)*cos(theta) + c*cos(phi)) * C.j +
        (a*cos(theta) - b*sin(theta)) * C.k)
    # vector expressed in two coordinate systems
    v5 = a * C.i + b * S.j
    assert express(v5, C) == (
        (a + b * cos(phi) * cos(theta)) * C.i +
        (b * sin(phi) * cos(theta)) * C.j +
        (-b * sin(theta)) * C.k)


def test_express_from_cartesian_to_spherical():
    C = CoordSys3D("C")
    S = C.create_new("S", transformation="spherical")
    r, theta, phi = S.base_scalars()
    a, b, c = symbols("a:c")

    # vector along the x-axis
    v1 = 1 * C.i
    assert express(v1, S) == (
        (sin(theta) * cos(phi)) * S.i +
        (cos(theta) * cos(phi)) * S.j +
        (-sin(phi)) * S.k)
    # vector along the y-axis
    v2 = 1 * C.j
    assert express(v2, S) == (
        (sin(theta) * sin(phi)) * S.i +
        (cos(theta) * sin(phi)) * S.j +
        (cos(phi)) * S.k)
    # vector along the z-axis
    v3 = 1 * C.k
    assert express(v3, S) == cos(theta) * S.i - sin(theta) * S.j
    # arbitrary vector expressed in cartesian coordinates
    v4 = a * C.i + b * C.j + c * C.k
    assert express(v4, S) == (
        (a*cos(phi)*sin(theta) + b*sin(phi)*sin(theta) + c*cos(theta)) * S.i +
        (a*cos(phi)*cos(theta) + b*sin(phi)*cos(theta) - c*sin(theta)) * S.j +
        (-a*sin(phi) + b*cos(phi)) * S.k)
    # vector expressed in two coordinate systems
    v5 = a * C.i + b * S.j
    assert express(v5, S) == (
        (a * cos(phi) * sin(theta)) * S.i +
        (a * cos(phi) * cos(theta) + b) * S.j +
        (-a * sin(phi)) * S.k)


def test_express_from_cylindrical_to_cartesian():
    Cart = CoordSys3D("Cart")
    Cyl = Cart.create_new("Cyl", transformation="cylindrical")
    r, theta, zeta = Cyl.base_scalars()
    a, b, c = symbols("a:c")

    # vector along the radial direction
    v1 = 1 * Cyl.i
    assert express(v1, Cart) == cos(theta) * Cart.i + sin(theta) * Cart.j
    # vector along the theta direction
    v2 = 1 * Cyl.j
    assert express(v2, Cart) == -sin(theta) * Cart.i + cos(theta) * Cart.j
    # vector along the z direction
    v3 = 1 * Cyl.k
    assert express(v3, Cart) == Cart.k
    # arbitrary vector expressed in cylindrical coordinates
    v4 = a * Cyl.i + b * Cyl.j + c * Cyl.k
    assert express(v4, Cart) == (
        (a*cos(theta) - b*sin(theta)) * Cart.i +
        (a*sin(theta) + b*cos(theta)) * Cart.j +
        c * Cart.k)
    # vector expressed in two coordinate systems
    v5 = a * Cart.i + b * Cyl.j
    assert express(v5, Cart) == (
        (a - b * sin(theta)) * Cart.i +
        (b * cos(theta)) * Cart.j)


def test_express_from_cartesian_to_cylindrical():
    Cart = CoordSys3D("Cart")
    Cyl = Cart.create_new("Cyl", transformation="cylindrical")
    r, theta, zeta = Cyl.base_scalars()
    a, b, c = symbols("a:c")

    # vector along the x-axis
    v1 = 1 * Cart.i
    assert express(v1, Cyl) == cos(theta) * Cyl.i - sin(theta) * Cyl.j
    # vector along the y-axis
    v2 = 1 * Cart.j
    assert express(v2, Cyl) == sin(theta) * Cyl.i + cos(theta) * Cyl.j
    # vector along the z-axis
    v3 = 1 * Cart.k
    assert express(v3, Cyl) == Cyl.k
    # arbitrary vector expressed in cartesian coordinates
    v4 = a * Cart.i + b * Cart.j + c * Cart.k
    assert express(v4, Cyl) == (
        (a*cos(theta) + b*sin(theta)) * Cyl.i +
        (-a*sin(theta) + b*cos(theta)) * Cyl.j +
        c * Cyl.k)
    # vector expressed in two coordinate systems
    v5 = a * Cart.i + b * Cyl.j
    assert express(v5, Cyl) == (
        (a * cos(theta)) * Cyl.i +
        (b - a * sin(theta)) * Cyl.j)


def test_express_from_spherical_to_cylindrical():
    Cart = CoordSys3D("Cart")
    S = Cart.create_new("S", transformation="spherical")
    C = Cart.create_new("C", transformation="cylindrical")
    r_s, theta_s, phi_s = S.base_scalars()
    r_c, theta_c, z_c = C.base_scalars()
    a, b, c = symbols("a, b, c")

    # The following tests uses simplify in order to show `theta_c - phis_s`,
    # which is a reminder to the user that the azimuthal angle
    # of S and C are the same, phi_s=theta_c, [0, 2*pi[

    # vector along the radial direction
    v1 = S.i
    r1 = express(v1, C)
    assert r1.simplify() == (
        (sin(theta_s) * cos(theta_c - phi_s)) * C.i +
        (-sin(theta_s) * sin(theta_c - phi_s)) * C.j +
        (cos(theta_s)) * C.k)
    assert r1.subs(phi_s, theta_c).simplify() == (
        sin(theta_s) * C.i + cos(theta_s) * C.k)

    # vector along the polar direction
    v2 = S.j
    r2 = express(v2, C)
    assert r2.simplify() == (
        (cos(theta_s) * cos(theta_c - phi_s)) * C.i +
        (-cos(theta_s) * sin(theta_c - phi_s)) * C.j +
        (-sin(theta_s)) * C.k)
    assert r2.subs(phi_s, theta_c).simplify() == (
        cos(theta_s) * C.i - sin(theta_s) * C.k)

    # vector along the azimuthal direction
    v3 = S.k
    r3 = express(v3, C)
    assert r3 == (
        (
            sin(C.theta)*cos(S.phi) -
            sin(S.phi)*cos(C.theta)
        ) * C.i +
        (
            sin(C.theta)*sin(S.phi) +
            cos(C.theta)*cos(S.phi)
        ) * C.j)
    assert r3.simplify() == (
        (sin(theta_c - phi_s)) * C.i +
        (cos(theta_c - phi_s)) * C.j)
    assert r3.subs(phi_s, theta_c).simplify() == C.j

    # vector along arbitrary direction
    v4 = a * S.i + b * S.j + c * S.k
    r4 = express(v4, C)
    assert r4.equals(
        (
            a*sin(S.theta)*cos(C.theta - S.phi) +
            b*cos(S.theta)*cos(C.theta - S.phi) +
            c*sin(C.theta - S.phi)
        ) * C.i +
        (
            -a*sin(S.theta)*sin(C.theta - S.phi) -
            b*sin(C.theta - S.phi)*cos(S.theta) +
            c*cos(C.theta - S.phi)
        ) * C.j +
        (
            a*cos(S.theta) - b*sin(S.theta)
        ) * C.k)
    assert r4.subs(phi_s, theta_c).equals(
        (a * sin(theta_s) + b * cos(theta_s)) * C.i +
        c * C.j +
        (a * cos(theta_s) - b * sin(theta_s)) * C.k)


def test_express_from_cylindrical_to_spherical():
    Cart = CoordSys3D("Cart")
    S = Cart.create_new("S", transformation="spherical")
    C = Cart.create_new("C", transformation="cylindrical")
    r_s, theta_s, phi_s = S.base_scalars()
    r_c, theta_c, z_c = C.base_scalars()
    a, b, c = symbols("a, b, c")

    # The following tests uses simplify in order to show `theta_c - phis_s`,
    # which is a reminder to the user that the azimuthal angle
    # of S and C are the same, phi_s=theta_c, [0, 2*pi[

    # vector along the radial direction
    v1 = C.i
    r1 = express(v1, S)
    assert r1.simplify() == (
        (sin(theta_s) * cos(theta_c - phi_s)) * S.i +
        (cos(theta_s) * cos(theta_c - phi_s)) * S.j +
        (sin(theta_c - phi_s)) * S.k)
    assert r1.subs(phi_s, theta_c).simplify() == (
        sin(theta_s) * S.i + cos(theta_s) * S.j)

    # vector along the azimuthal direction
    v2 = C.j
    r2 = express(v2, S)
    assert r2.simplify() == (
        (-sin(theta_s) * sin(theta_c - phi_s)) * S.i +
        (-cos(theta_s) * sin(theta_c - phi_s)) * S.j +
        (cos(theta_c - phi_s)) * S.k
    )
    assert r2.subs(phi_s, theta_c).simplify() == S.k

    # vector along the z-axis
    v3 = C.k
    r3 = express(v3, S)
    assert r3 == cos(theta_s) * S.i - sin(theta_s) * S.j
    assert r3.subs(phi_s, theta_c) == cos(theta_s) * S.i - sin(theta_s) * S.j

    # vector along arbitrary direction
    v4 = a * C.i + b * C.j + c * C.k
    r4 = express(v4, S)
    assert r4.simplify() == (
        (
            a*sin(S.theta)*cos(C.theta - S.phi) -
            b*sin(S.theta)*sin(C.theta - S.phi) +
            c*cos(S.theta)
        ) * S.i +
        (
            a*cos(S.theta)*cos(C.theta - S.phi) -
            b*sin(C.theta - S.phi)*cos(S.theta) -
            c*sin(S.theta)
        ) * S.j +
        (
            a*sin(C.theta - S.phi) + b*cos(C.theta - S.phi)
        ) * S.k)
    assert r4.subs(phi_s, theta_c).simplify() == (
        (a * sin(theta_s) + c * cos(theta_s)) * S.i +
        (a * cos(theta_s) - c * sin(theta_s)) * S.j +
        b * S.k)


def test_express_scaled_and_reflected_to_cartesian():
    C = CoordSys3D("C")
    A = C.create_new("A", transformation=lambda x, y, z: (y, 2*x, z))
    a, b, c = symbols("a, b, c")

    assert express(A.i, C) == 2 * C.j
    assert express(A.j, C) == C.i
    assert express(A.k, C) == C.k
    assert express(a * A.i + b * A.j + c * A.k, C) == (
        b * C.i + 2 * a * C.j + c * C.k)


def test_express_cartesian_to_scaled_and_reflected():
    C = CoordSys3D("C")
    A = C.create_new("A", transformation=lambda x, y, z: (y, 2*x, z))
    a, b, c = symbols("a, b, c")

    assert express(C.i, A) == A.j
    assert express(C.j, A) == A.i / 2
    assert express(C.k, A) == A.k
    assert express(a * C.i + b * C.j + c * C.k, A) == (
        b / 2 * A.i + a * A.j + c * A.k)


def test_cartesian_flipped_to_cartesian():
    C1 = CoordSys3D("C1")
    C2 = C1.create_new("C2", transformation=lambda x, y, z: (y, x, z))
    a, b, c = symbols("a, b, c")

    assert express(C2.i, C1) == C1.j
    assert express(C2.j, C1) == C1.i
    assert express(C2.k, C1) == C1.k
    assert express(a * C2.i + b * C2.j + c * C2.k, C1) == (
        b * C1.i + a * C1.j + c * C1.k)


def test_cartesian_to_cartesian_flipped():
    C1 = CoordSys3D("C1")
    C2 = C1.create_new("C2", transformation=lambda x, y, z: (y, x, z))
    a, b, c = symbols("a, b, c")

    assert express(C1.i, C2) == C2.j
    assert express(C1.j, C2) == C2.i
    assert express(C1.k, C2) == C2.k
    assert express(a * C1.i + b * C1.j + c * C1.k, C2) == (
        b * C2.i + a * C2.j + c * C2.k)

def test_two_cartesian_flipped_systems():
    A = CoordSys3D("A", transformation=lambda x, y, z: (y, x, z))
    B = A.create_new("B", transformation=lambda x, y, z: (z, y, x))
    a, b, c = symbols("a:c")

    assert express(a * B.i + b * B.j + c * B.k, A) == (
        c * A.i + b * A.j + a * A.k)
    assert express(a * A.i + b * A.j + c * A.k, B) == (
        c * B.i + b * B.j + a * B.k)

    C = CoordSys3D("C")
    A = C.create_new("A", transformation=lambda x, y, z: (y, x, z))
    B = C.create_new("B", transformation=lambda x, y, z: (z, y, x))
    assert express(a * B.i + b * B.j + c * B.k, A) == (
        b * A.i + c * A.j + a * A.k)
    assert express(a * A.i + b * A.j + c * A.k, B) == (
        c * B.i + a * B.j + b * B.k)


def test_express_path_of_systems():
    # C1 -> C2 (translate) -> C3 (rotate about k by alpha) -> S (spherical)
    # Note that the transformation matrix doesn't consider translation,
    # only rotation and scaling
    p, alpha = symbols("p, alpha")
    C1 = CoordSys3D("C1")
    C2 = C1.orient_new_axis("C2", alpha, C1.k)
    C3 = C2.locate_new("C3", p * C2.i)
    S = C3.create_new("S", transformation="spherical")
    r, theta, phi = S.base_scalars()

    # vector along the radial direction
    v1 = S.i
    assert express(v1, C1) == (
        (
            -sin(phi) * sin(theta) * sin(alpha) +
            cos(phi) * sin(theta) * cos(alpha)
        ) * C1.i +
        (
            sin(phi) * sin(theta) * cos(alpha) +
            cos(phi) * sin(theta) * sin(alpha)
        ) * C1.j +
        cos(theta) * C1.k)
    assert express(v1, C3) == (
        (sin(theta) * cos(phi)) * C3.i +
        (sin(theta) * sin(phi)) * C3.j +
        cos(theta) * C3.k)

    # vector along the polar direction
    v2 = 1 * S.j
    assert express(v2, C1) == (
        (
            -sin(phi) * cos(theta) * sin(alpha) +
            cos(phi) * cos(theta) * cos(alpha)
        ) * C1.i +
        (
            sin(phi) * cos(theta) * cos(alpha) +
            cos(phi) * cos(theta) * sin(alpha)
        ) * C1.j +
        (-sin(theta)) * C1.k)
    assert express(v2, C3) == (
        (cos(phi) * cos(theta)) * C3.i +
        (sin(phi) * cos(theta)) * C3.j +
        (-sin(theta)) * C3.k)

    # vector along the azimuthal direction
    v3 = 1 * S.k
    assert express(v3, C1) == (
        (-sin(phi) * cos(alpha) - sin(alpha) * cos(phi)) * C1.i +
        (-sin(phi) * sin(alpha) + cos(phi) * cos(alpha)) * C1.j)
    assert express(v3, C3) == (
        -sin(phi) * C3.i +
        cos(phi) * C3.j)


def test_express_with_variables():
    A = CoordSys3D("A")
    C = A.create_new('C', transformation='cylindrical')

    assert express(A.x, C, variables=False) == A.x
    assert express(A.x, C, variables=True) == C.r * cos(C.theta)
    assert express(A.x*A.j, C, variables=False) == (
        A.x * sin(C.theta) * C.i + A.x * cos(C.theta) * C.j)
    res = express(A.x*A.j, C, variables=True)
    assert res == (
        C.r * sin(C.theta) * cos(C.theta) * C.i +
        C.r * cos(C.theta)**2 * C.j)
    assert express(res, A, variables=False) == (
        (C.r * sin(C.theta)**2 * cos(C.theta) + C.r * cos(C.theta)**3) * A.j)
    assert express(res, A, variables=True) == (
        (A.x**3 / (A.x**2 + A.y**2) + A.x * A.y**2 / (A.x**2 + A.y**2)) * A.j)

    assert express(C.r, A, variables=False) == C.r
    assert express(C.r, A, variables=True) == sqrt(A.x**2 + A.y**2)
    assert express(C.r * C.i, A, variables=False) == (
        C.r * cos(C.theta) * A.i + C.r * sin(C.theta) * A.j)
    res = express(C.r * C.i, A, variables=True)
    assert res == A.x * A.i + A.y * A.j
    assert express(res, C, variables=False) == (
        (A.x * cos(C.theta) + A.y * sin(C.theta)) * C.i +
        (-A.x * sin(C.theta) + A.y * cos(C.theta)) * C.j)
    assert express(res, C, variables=True) == (
        (C.r * sin(C.theta)**2 + C.r * cos(C.theta)**2) * C.i)
    assert express(res, C, variables=True).simplify() == C.r * C.i
