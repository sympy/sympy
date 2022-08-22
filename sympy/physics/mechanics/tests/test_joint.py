from sympy.core.function import expand_mul
from sympy.core.numbers import pi
from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.matrices.dense import Matrix
from sympy.core.backend import _simplify_matrix
from sympy.core.symbol import symbols
from sympy.physics.mechanics import dynamicsymbols, Body, PinJoint, PrismaticJoint
from sympy.physics.mechanics.joint import Joint
from sympy.physics.vector import Vector, ReferenceFrame, Point
from sympy.testing.pytest import (raises, ignore_warnings, XFAIL,
                                  SymPyDeprecationWarning,
                                  warns_deprecated_sympy)

t = dynamicsymbols._t # type: ignore


def _generate_body(interframe=False):
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    P = Body('P', frame=N)
    C = Body('C', frame=A)
    if interframe:
        Pint, Cint = ReferenceFrame('P_int'), ReferenceFrame('C_int')
        Pint.orient_axis(N, N.x, pi)
        Cint.orient_axis(A, A.y, -pi / 2)
        return N, A, P, C, Pint, Cint
    return N, A, P, C


def test_Joint():
    parent = Body('parent')
    child = Body('child')
    raises(TypeError, lambda: Joint('J', parent, child))


def test_pinjoint():
    P = Body('P')
    C = Body('C')
    l, m = symbols('l m')
    theta, omega = dynamicsymbols('theta_J, omega_J')
    Pj = PinJoint('J', P, C)
    assert Pj.name == 'J'
    assert Pj.parent == P
    assert Pj.child == C
    assert Pj.coordinates == [theta]
    assert Pj.speeds == [omega]
    assert Pj.kdes == [omega - theta.diff(t)]
    assert Pj.joint_axis == P.frame.x
    assert Pj.child_point.pos_from(C.masscenter) == Vector(0)
    assert Pj.parent_point.pos_from(P.masscenter) == Vector(0)
    assert Pj.parent_point.pos_from(Pj._child_point) == Vector(0)
    assert C.masscenter.pos_from(P.masscenter) == Vector(0)
    assert Pj.parent_interframe == P.frame
    assert Pj.child_interframe == C.frame
    assert Pj.__str__() == 'PinJoint: J  parent: P  child: C'

    P1 = Body('P1')
    C1 = Body('C1')
    Pint = ReferenceFrame('P_int')
    Pint.orient_axis(P1.frame, P1.y, pi / 2)
    J1 = PinJoint('J1', P1, C1, parent_point=l*P1.frame.x,
                  child_point=m*C1.frame.y, joint_axis=P1.frame.z,
                  parent_interframe=Pint)
    assert J1._joint_axis == P1.frame.z
    assert J1._child_point.pos_from(C1.masscenter) == m * C1.frame.y
    assert J1._parent_point.pos_from(P1.masscenter) == l * P1.frame.x
    assert J1._parent_point.pos_from(J1._child_point) == Vector(0)
    assert (P1.masscenter.pos_from(C1.masscenter) ==
            -l*P1.frame.x + m*C1.frame.y)
    assert J1.parent_interframe == Pint
    assert J1.child_interframe == C1.frame

    q, u = dynamicsymbols('q, u')
    N, A, P, C, Pint, Cint = _generate_body(True)
    parent_point = P.masscenter.locatenew('parent_point', N.x + N.y)
    child_point = C.masscenter.locatenew('child_point', C.y + C.z)
    J = PinJoint('J', P, C, q, u, parent_point=parent_point,
                 child_point=child_point, parent_interframe=Pint,
                 child_interframe=Cint, joint_axis=N.z)
    assert J.joint_axis == N.z
    assert J.parent_point.vel(N) == 0
    assert J.parent_point == parent_point
    assert J.child_point == child_point
    assert J.child_point.pos_from(P.masscenter) == N.x + N.y
    assert J.parent_point.pos_from(C.masscenter) == C.y + C.z
    assert C.masscenter.pos_from(P.masscenter) == N.x + N.y - C.y - C.z
    assert C.masscenter.vel(N).express(N) == (u * sin(q) - u * cos(q)) * N.x + (
            -u * sin(q) - u * cos(q)) * N.y
    assert J.parent_interframe == Pint
    assert J.child_interframe == Cint


def test_pin_joint_double_pendulum():
    q1, q2 = dynamicsymbols('q1 q2')
    u1, u2 = dynamicsymbols('u1 u2')
    m, l = symbols('m l')
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    B = ReferenceFrame('B')
    C = Body('C', frame=N)  # ceiling
    PartP = Body('P', frame=A, mass=m)
    PartR = Body('R', frame=B, mass=m)

    J1 = PinJoint('J1', C, PartP, speeds=u1, coordinates=q1,
                  child_point=-l*A.x, joint_axis=C.frame.z)
    J2 = PinJoint('J2', PartP, PartR, speeds=u2, coordinates=q2,
                  child_point=-l*B.x, joint_axis=PartP.frame.z)

    # Check orientation
    assert N.dcm(A) == Matrix([[cos(q1), -sin(q1), 0],
                               [sin(q1), cos(q1), 0], [0, 0, 1]])
    assert A.dcm(B) == Matrix([[cos(q2), -sin(q2), 0],
                               [sin(q2), cos(q2), 0], [0, 0, 1]])
    assert _simplify_matrix(N.dcm(B)) == Matrix([[cos(q1 + q2), -sin(q1 + q2), 0],
                                                 [sin(q1 + q2), cos(q1 + q2), 0],
                                                 [0, 0, 1]])

    # Check Angular Velocity
    assert A.ang_vel_in(N) == u1 * N.z
    assert B.ang_vel_in(A) == u2 * A.z
    assert B.ang_vel_in(N) == u1 * N.z + u2 * A.z

    # Check kde
    assert J1.kdes == [u1 - q1.diff(t)]
    assert J2.kdes == [u2 - q2.diff(t)]

    # Check Linear Velocity
    assert PartP.masscenter.vel(N) == l*u1*A.y
    assert PartR.masscenter.vel(A) == l*u2*B.y
    assert PartR.masscenter.vel(N) == l*u1*A.y + l*(u1 + u2)*B.y


def test_pin_joint_chaos_pendulum():
    mA, mB, lA, lB, h = symbols('mA, mB, lA, lB, h')
    theta, phi, omega, alpha = dynamicsymbols('theta phi omega alpha')
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    B = ReferenceFrame('B')
    lA = (lB - h / 2) / 2
    lC = (lB/2 + h/4)
    rod = Body('rod', frame=A, mass=mA)
    plate = Body('plate', mass=mB, frame=B)
    C = Body('C', frame=N)
    J1 = PinJoint('J1', C, rod, coordinates=theta, speeds=omega,
                  child_point=lA*A.z, joint_axis=N.y)
    J2 = PinJoint('J2', rod, plate, coordinates=phi, speeds=alpha,
                  parent_point=lC*A.z, joint_axis=A.z)

    # Check orientation
    assert A.dcm(N) == Matrix([[cos(theta), 0, -sin(theta)],
                               [0, 1, 0],
                               [sin(theta), 0, cos(theta)]])
    assert A.dcm(B) == Matrix([[cos(phi), -sin(phi), 0],
                               [sin(phi), cos(phi), 0],
                               [0, 0, 1]])
    assert B.dcm(N) == Matrix([
        [cos(phi)*cos(theta), sin(phi), -sin(theta)*cos(phi)],
        [-sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta)],
        [sin(theta), 0, cos(theta)]])

    # Check Angular Velocity
    assert A.ang_vel_in(N) == omega*N.y
    assert A.ang_vel_in(B) == -alpha*A.z
    assert N.ang_vel_in(B) == -omega*N.y - alpha*A.z

    # Check kde
    assert J1.kdes == [omega - theta.diff(t)]
    assert J2.kdes == [alpha - phi.diff(t)]

    # Check pos of masscenters
    assert C.masscenter.pos_from(rod.masscenter) == lA*A.z
    assert rod.masscenter.pos_from(plate.masscenter) == - lC * A.z

    # Check Linear Velocities
    assert rod.masscenter.vel(N) == (h/4 - lB/2)*omega*A.x
    assert plate.masscenter.vel(N) == ((h/4 - lB/2)*omega +
                                       (h/4 + lB/2)*omega)*A.x


@XFAIL
def test_pinjoint_interframe():
    q, u = dynamicsymbols('q, u')
    # Check not connected
    N, A, P, C = _generate_body()
    Pint, Cint = ReferenceFrame('Pint'), ReferenceFrame('Cint')
    raises(ValueError, lambda: PinJoint('J', P, C, parent_interframe=Pint))
    raises(ValueError, lambda: PinJoint('J', P, C, child_interframe=Cint))
    # Check not fixed interframe
    Pint.orient_axis(N, N.z, q)
    Cint.orient_axis(A, A.z, q)
    raises(ValueError, lambda: PinJoint('J', P, C, parent_interframe=Pint))
    raises(ValueError, lambda: PinJoint('J', P, C, child_interframe=Cint))
    # Check only parent_interframe
    N, A, P, C = _generate_body()
    Pint = ReferenceFrame('Pint')
    Pint.orient_body_fixed(N, (pi / 4, pi, pi / 3), 'xyz')
    PinJoint('J', P, C, q, u, parent_point=N.x, child_point=-C.y,
             parent_interframe=Pint, joint_axis=C.x)
    assert _simplify_matrix(N.dcm(A)) == Matrix([
        [-1 / 2, sqrt(3) * cos(q) / 2, -sqrt(3) * sin(q) / 2],
        [sqrt(6) / 4, sqrt(2) * (2 * sin(q) + cos(q)) / 4,
         sqrt(2) * (-sin(q) + 2 * cos(q)) / 4],
        [sqrt(6) / 4, sqrt(2) * (-2 * sin(q) + cos(q)) / 4,
         -sqrt(2) * (sin(q) + 2 * cos(q)) / 4]])
    assert A.ang_vel_in(N) == u * Pint.x
    assert C.masscenter.pos_from(P.masscenter) == N.x + A.y
    assert C.masscenter.vel(N) == u * A.z
    # Check only child_interframe
    N, A, P, C = _generate_body()
    Cint = ReferenceFrame('Cint')
    Cint.orient_body_fixed(A, (2 * pi / 3, -pi, pi / 2), 'xyz')
    PinJoint('J', P, C, q, u, parent_point=-N.z, child_point=C.x,
             child_interframe=Cint, joint_axis=P.x + P.z)
    assert _simplify_matrix(N.dcm(A)) == Matrix([
        [-sqrt(2) * sin(q) / 2,
         -sqrt(3) * (cos(q) - 1) / 4 - cos(q) / 4 - S(1) / 4,
         sqrt(3) * (cos(q) + 1) / 4 - cos(q) / 4 + S(1) / 4],
        [cos(q), (sqrt(2) + sqrt(6)) * -sin(q) / 4,
         (-sqrt(2) + sqrt(6)) * sin(q) / 4],
        [sqrt(2) * sin(q) / 2,
         sqrt(3) * (cos(q) + 1) / 4 + cos(q) / 4 - S(1) / 4,
         sqrt(3) * (1 - cos(q)) / 4 + cos(q) / 4 + S(1) / 4]])
    assert A.ang_vel_in(N) == sqrt(2) * u / 2 * N.x + sqrt(2) * u / 2 * N.z
    assert C.masscenter.pos_from(P.masscenter) == - N.z - A.x
    assert C.masscenter.vel(N).simplify() == (
        -sqrt(6) - sqrt(2)) * u / 4 * A.y + (
               -sqrt(2) + sqrt(6)) * u / 4 * A.z
    # Check combination
    N, A, P, C = _generate_body()
    Pint, Cint = ReferenceFrame('Pint'), ReferenceFrame('Cint')
    Pint.orient_body_fixed(N, (-pi / 2, pi, pi / 2), 'xyz')
    Cint.orient_body_fixed(A, (2 * pi / 3, -pi, pi / 2), 'xyz')
    PinJoint('J', P, C, q, u, parent_point=N.x - N.y, child_point=-C.z,
             parent_interframe=Pint, child_interframe=Cint,
             joint_axis=Cint.x + Cint.z)
    assert _simplify_matrix(N.dcm(A)) == Matrix([
        [cos(q), (sqrt(2) + sqrt(6)) * -sin(q) / 4,
         (-sqrt(2) + sqrt(6)) * sin(q) / 4],
        [-sqrt(2) * sin(q) / 2,
         -sqrt(3) * (cos(q) + 1) / 4 - cos(q) / 4 + S(1) / 4,
         sqrt(3) * (cos(q) - 1) / 4 - cos(q) / 4 - S(1) / 4],
        [sqrt(2) * sin(q) / 2,
         sqrt(3) * (cos(q) - 1) / 4 + cos(q) / 4 + S(1) / 4,
         -sqrt(3) * (cos(q) + 1) / 4 + cos(q) / 4 - S(1) / 4]])
    assert A.ang_vel_in(N) == sqrt(2) * u / 2 * Pint.x + sqrt(
        2) * u / 2 * Pint.z
    assert C.masscenter.pos_from(P.masscenter) == N.x - N.y + A.z
    assert C.masscenter.vel(N).simplify() == (-sqrt(2) + sqrt(6)) * u / 4 * A.x


def test_pinjoint_joint_axis():
    q, u = dynamicsymbols('q, u')
    # Check parent as reference
    N, A, P, C, Pint, Cint = _generate_body(True)
    pin = PinJoint('J', P, C, q, u, parent_interframe=Pint,
                   child_interframe=Cint, joint_axis=P.y)
    assert pin.joint_axis == P.y
    assert N.dcm(A) == Matrix([[sin(q), 0, cos(q)], [0, -1, 0],
                               [cos(q), 0, -sin(q)]])
    # Check parent_interframe as reference
    N, A, P, C, Pint, Cint = _generate_body(True)
    pin = PinJoint('J', P, C, q, u, parent_interframe=Pint,
                   child_interframe=Cint, joint_axis=Pint.y)
    assert pin.joint_axis == Pint.y
    assert N.dcm(A) == Matrix([[-sin(q), 0, cos(q)], [0, -1, 0],
                               [cos(q), 0, sin(q)]])
    # Check child_interframe as reference
    N, A, P, C, Pint, Cint = _generate_body(True)
    pin = PinJoint('J', P, C, q, u, parent_interframe=Pint,
                   child_interframe=Cint, joint_axis=Cint.y)
    assert pin.joint_axis == Cint.y
    assert N.dcm(A) == Matrix([[-sin(q), 0, cos(q)], [0, -1, 0],
                               [cos(q), 0, sin(q)]])
    # Check child as reference
    N, A, P, C, Pint, Cint = _generate_body(True)
    pin = PinJoint('J', P, C, q, u, parent_interframe=Pint,
                   child_interframe=Cint, joint_axis=C.y)
    assert pin.joint_axis == C.y
    assert N.dcm(A) == Matrix([[-sin(q), 0, cos(q)], [0, -1, 0],
                               [cos(q), 0, sin(q)]])
    # Check time varying axis
    N, A, P, C, Pint, Cint = _generate_body(True)
    raises(ValueError, lambda: PinJoint('J', P, C,
                                        joint_axis=cos(q) * N.x + sin(q) * N.y))
    # Check some invalid combinations
    raises(ValueError, lambda: PinJoint('J', P, C, joint_axis=P.x + C.y))
    raises(ValueError, lambda: PinJoint(
        'J', P, C, parent_interframe=Pint, child_interframe=Cint,
        joint_axis=Pint.x + C.y))
    raises(ValueError, lambda: PinJoint(
        'J', P, C, parent_interframe=Pint, child_interframe=Cint,
        joint_axis=P.x + Cint.y))
    # Check valid special combination
    N, A, P, C, Pint, Cint = _generate_body(True)
    PinJoint('J', P, C, parent_interframe=Pint, child_interframe=Cint,
             joint_axis=Pint.x + P.y)
    N, A, P, C, Pint, Cint = _generate_body(True)
    PinJoint('J', P, C, parent_interframe=Pint, child_interframe=Cint,
             joint_axis=Cint.x + C.y)
    # Check invalid zero vector
    raises(Exception, lambda: PinJoint(
        'J', P, C, parent_interframe=Pint, child_interframe=Cint,
        joint_axis=Vector(0)))
    raises(ValueError, lambda: PinJoint(
        'J', P, C, parent_interframe=Pint, child_interframe=Cint,
        joint_axis=C.z - Cint.x))


def test_pinjoint_arbitrary_axis():
    theta, omega = dynamicsymbols('theta_J, omega_J')

    # When the bodies are attached though masscenters but axess are opposite.
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PinJoint('J', P, C, child_axis=-A.x)

    assert (-A.x).angle_between(N.x) == 0
    assert -A.x.express(N) == N.x
    assert A.dcm(N) == Matrix([[-1, 0, 0],
                            [0, -cos(theta), -sin(theta)],
                            [0, -sin(theta), cos(theta)]])
    assert A.ang_vel_in(N) == omega*N.x
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    assert C.masscenter.pos_from(P.masscenter) == 0
    assert C.masscenter.pos_from(P.masscenter).express(N).simplify() == 0
    assert C.masscenter.vel(N) == 0

    # When axes are different and parent joint is at masscenter but child joint
    # is at a unit vector from child masscenter.
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PinJoint('J', P, C, child_axis=A.y, child_point=A.x)

    assert A.y.angle_between(N.x) == 0  # Axis are aligned
    assert A.y.express(N) == N.x
    assert A.dcm(N) == Matrix([[0, -cos(theta), -sin(theta)],
                               [1, 0, 0],
                               [0, -sin(theta), cos(theta)]])
    assert A.ang_vel_in(N) == omega*N.x
    assert A.ang_vel_in(N).express(A) == omega * A.y
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    assert A.ang_vel_in(N).cross(A.y) == 0
    assert C.masscenter.vel(N) == omega*A.z
    assert C.masscenter.pos_from(P.masscenter) == -A.x
    assert (C.masscenter.pos_from(P.masscenter).express(N).simplify() ==
            cos(theta)*N.y + sin(theta)*N.z)
    assert C.masscenter.vel(N).angle_between(A.x) == pi/2

    # Similar to previous case but wrt parent body
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PinJoint('J', P, C, parent_axis=N.y, parent_point=N.x)

    assert N.y.angle_between(A.x) == 0  # Axis are aligned
    assert N.y.express(A) == A.x
    assert A.dcm(N) == Matrix([[0, 1, 0],
                               [-cos(theta), 0, sin(theta)],
                               [sin(theta), 0, cos(theta)]])
    assert A.ang_vel_in(N) == omega*N.y
    assert A.ang_vel_in(N).express(A) == omega*A.x
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    angle = A.ang_vel_in(N).angle_between(A.x)
    assert angle.xreplace({omega: 1}) == 0
    assert C.masscenter.vel(N) == 0
    assert C.masscenter.pos_from(P.masscenter) == N.x

    # Both joint pos id defined but different axes
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PinJoint('J', P, C, parent_point=N.x, child_point=A.x,
                 child_axis=A.x+A.y)
    assert expand_mul(N.x.angle_between(A.x + A.y)) == 0  # Axis are aligned
    assert (A.x + A.y).express(N).simplify() == sqrt(2)*N.x
    assert _simplify_matrix(A.dcm(N)) == Matrix([
        [sqrt(2)/2, -sqrt(2)*cos(theta)/2, -sqrt(2)*sin(theta)/2],
        [sqrt(2)/2, sqrt(2)*cos(theta)/2, sqrt(2)*sin(theta)/2],
        [0, -sin(theta), cos(theta)]])
    assert A.ang_vel_in(N) == omega*N.x
    assert (A.ang_vel_in(N).express(A).simplify() ==
            (omega*A.x + omega*A.y)/sqrt(2))
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    angle = A.ang_vel_in(N).angle_between(A.x + A.y)
    assert angle.xreplace({omega: 1}) == 0
    assert C.masscenter.vel(N).simplify() == (omega * A.z)/sqrt(2)
    assert C.masscenter.pos_from(P.masscenter) == N.x - A.x
    assert (C.masscenter.pos_from(P.masscenter).express(N).simplify() ==
            (1 - sqrt(2)/2)*N.x + sqrt(2)*cos(theta)/2*N.y +
            sqrt(2)*sin(theta)/2*N.z)
    assert (C.masscenter.vel(N).express(N).simplify() ==
            -sqrt(2)*omega*sin(theta)/2*N.y + sqrt(2)*omega*cos(theta)/2*N.z)
    assert C.masscenter.vel(N).angle_between(A.x) == pi/2

    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PinJoint('J', P, C, parent_point=N.x, child_point=A.x,
                 child_axis=A.x+A.y-A.z)
    assert expand_mul(N.x.angle_between(A.x + A.y - A.z)) == 0  # Axis aligned
    assert (A.x + A.y - A.z).express(N).simplify() == sqrt(3)*N.x
    assert _simplify_matrix(A.dcm(N)) == Matrix([
        [sqrt(3)/3, -sqrt(6)*sin(theta + pi/4)/3,
         sqrt(6)*cos(theta + pi/4)/3],
        [sqrt(3)/3, sqrt(6)*cos(theta + pi/12)/3,
         sqrt(6)*sin(theta + pi/12)/3],
        [-sqrt(3)/3, sqrt(6)*cos(theta + 5*pi/12)/3,
         sqrt(6)*sin(theta + 5*pi/12)/3]])
    assert A.ang_vel_in(N) == omega*N.x
    assert A.ang_vel_in(N).express(A).simplify() == (omega*A.x + omega*A.y -
                                                     omega*A.z)/sqrt(3)
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    angle = A.ang_vel_in(N).angle_between(A.x + A.y-A.z)
    assert angle.xreplace({omega: 1}) == 0
    assert C.masscenter.vel(N).simplify() == (omega*A.y + omega*A.z)/sqrt(3)
    assert C.masscenter.pos_from(P.masscenter) == N.x - A.x
    assert (C.masscenter.pos_from(P.masscenter).express(N).simplify() ==
            (1 - sqrt(3)/3)*N.x + sqrt(6)*sin(theta + pi/4)/3*N.y -
            sqrt(6)*cos(theta + pi/4)/3*N.z)
    assert (C.masscenter.vel(N).express(N).simplify() ==
            sqrt(6)*omega*cos(theta + pi/4)/3*N.y +
            sqrt(6)*omega*sin(theta + pi/4)/3*N.z)
    assert C.masscenter.vel(N).angle_between(A.x) == pi/2

    N, A, P, C = _generate_body()
    m, n = symbols('m n')
    with warns_deprecated_sympy():
        PinJoint('J', P, C, parent_point=m*N.x, child_point=n*A.x,
                 child_axis=A.x+A.y-A.z, parent_axis=N.x-N.y+N.z)
    angle = (N.x-N.y+N.z).angle_between(A.x+A.y-A.z)
    assert expand_mul(angle) == 0  # Axis are aligned
    assert ((A.x-A.y+A.z).express(N).simplify() ==
            (-4*cos(theta)/3 - S(1)/3)*N.x + (S(1)/3 - 4*sin(theta + pi/6)/3)*N.y +
            (4*cos(theta + pi/3)/3 - S(1)/3)*N.z)
    assert _simplify_matrix(A.dcm(N)) == Matrix([
        [S(1)/3 - 2*cos(theta)/3, -2*sin(theta + pi/6)/3 - S(1)/3,
         2*cos(theta + pi/3)/3 + S(1)/3],
        [2*cos(theta + pi/3)/3 + S(1)/3, 2*cos(theta)/3 - S(1)/3,
         2*sin(theta + pi/6)/3 + S(1)/3],
        [-2*sin(theta + pi/6)/3 - S(1)/3, 2*cos(theta + pi/3)/3 + S(1)/3,
         2*cos(theta)/3 - S(1)/3]])
    assert A.ang_vel_in(N) == (omega*N.x - omega*N.y + omega*N.z)/sqrt(3)
    assert A.ang_vel_in(N).express(A).simplify() == (omega*A.x + omega*A.y -
                                                     omega*A.z)/sqrt(3)
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    angle = A.ang_vel_in(N).angle_between(A.x+A.y-A.z)
    assert angle.xreplace({omega: 1}) == 0
    assert (C.masscenter.vel(N).simplify() ==
            sqrt(3)*n*omega/3*A.y + sqrt(3)*n*omega/3*A.z)
    assert C.masscenter.pos_from(P.masscenter) == m*N.x - n*A.x
    assert (C.masscenter.pos_from(P.masscenter).express(N).simplify() ==
            (m + n*(2*cos(theta) - 1)/3)*N.x + n*(2*sin(theta + pi/6) +
            1)/3*N.y - n*(2*cos(theta + pi/3) + 1)/3*N.z)
    assert (C.masscenter.vel(N).express(N).simplify() ==
            - 2*n*omega*sin(theta)/3*N.x + 2*n*omega*cos(theta + pi/6)/3*N.y +
            2*n*omega*sin(theta + pi/3)/3*N.z)
    assert C.masscenter.vel(N).dot(N.x - N.y + N.z).simplify() == 0


def test_pinjoint_pi():
    _, _, P, C = _generate_body()
    with ignore_warnings(SymPyDeprecationWarning):
        J = PinJoint('J', P, C, child_axis=-C.frame.x)
        assert J._generate_vector() == P.frame.z

        _, _, P, C = _generate_body()
        J = PinJoint('J', P, C, parent_axis=P.frame.y, child_axis=-C.frame.y)
        assert J._generate_vector() == P.frame.x

        _, _, P, C = _generate_body()
        J = PinJoint('J', P, C, parent_axis=P.frame.z, child_axis=-C.frame.z)
        assert J._generate_vector() == P.frame.y

        _, _, P, C = _generate_body()
        J = PinJoint('J', P, C, parent_axis=P.frame.x+P.frame.y,
                     child_axis=-C.frame.y-C.frame.x)
        assert J._generate_vector() == P.frame.z

        _, _, P, C = _generate_body()
        J = PinJoint('J', P, C, parent_axis=P.frame.y+P.frame.z,
                     child_axis=-C.frame.y-C.frame.z)
        assert J._generate_vector() == P.frame.x

        _, _, P, C = _generate_body()
        J = PinJoint('J', P, C, parent_axis=P.frame.x+P.frame.z,
                     child_axis=-C.frame.z-C.frame.x)
        assert J._generate_vector() == P.frame.y

        _, _, P, C = _generate_body()
        J = PinJoint('J', P, C, parent_axis=P.frame.x+P.frame.y+P.frame.z,
                     child_axis=-C.frame.x-C.frame.y-C.frame.z)
        assert J._generate_vector() == P.frame.y - P.frame.z


def test_pinjoint_axis():
    q, u = dynamicsymbols('q u')
    # Test default joint axis
    N, A, P, C, Pint, Cint = _generate_body(True)
    J = PinJoint('J', P, C, q, u, parent_interframe=Pint, child_interframe=Cint)
    assert J.joint_axis == Pint.x
    # Test for the same joint axis expressed in different frames
    N_R_A = Matrix([[0, sin(q), cos(q)],
                    [0, -cos(q), sin(q)],
                    [1, 0, 0]])
    N, A, P, C, Pint, Cint = _generate_body(True)
    PinJoint('J', P, C, q, u, parent_interframe=Pint, child_interframe=Cint,
             joint_axis=N.z)
    assert N.dcm(A) == N_R_A
    N, A, P, C, Pint, Cint = _generate_body(True)
    PinJoint('J', P, C, q, u, parent_interframe=Pint, child_interframe=Cint,
             joint_axis=-Pint.z)
    assert N.dcm(A) == N_R_A
    N, A, P, C, Pint, Cint = _generate_body(True)
    PinJoint('J', P, C, q, u, parent_interframe=Pint, child_interframe=Cint,
             joint_axis=-Cint.z)
    assert N.dcm(A) == N_R_A
    N, A, P, C, Pint, Cint = _generate_body(True)
    PinJoint('J', P, C, q, u, parent_interframe=Pint, child_interframe=Cint,
             joint_axis=A.x)
    assert N.dcm(A) == N_R_A
    # Test time varying joint axis
    N, A, P, C, Pint, Cint = _generate_body(True)
    raises(ValueError, lambda: PinJoint('J', P, C, joint_axis=q * N.z))


def test_locate_joint_pos():
    # Test Vector and default
    N, A, P, C = _generate_body()
    joint = PinJoint('J', P, C, parent_point=N.y + N.z)
    assert joint.parent_point.name == 'J_P_joint'
    assert joint.parent_point.pos_from(P.masscenter) == N.y + N.z
    assert joint.child_point == C.masscenter
    # Test Point objects
    N, A, P, C = _generate_body()
    parent_point = P.masscenter.locatenew('p', N.y + N.z)
    joint = PinJoint('J', P, C, parent_point=parent_point,
                     child_point=C.masscenter)
    assert joint.parent_point == parent_point
    assert joint.child_point == C.masscenter
    # Check invalid type
    N, A, P, C = _generate_body()
    raises(TypeError,
           lambda: PinJoint('J', P, C, parent_point=N.x.to_matrix(N)))
    # Test time varying positions
    q = dynamicsymbols('q')
    N, A, P, C = _generate_body()
    raises(ValueError, lambda: PinJoint('J', P, C, parent_point=q * N.x))
    N, A, P, C = _generate_body()
    child_point = C.masscenter.locatenew('p', q * A.y)
    raises(ValueError, lambda: PinJoint('J', P, C, child_point=child_point))
    # Test undefined position
    child_point = Point('p')
    raises(ValueError, lambda: PinJoint('J', P, C, child_point=child_point))


def test_locate_joint_frame():
    # Test rotated frame and default
    N, A, P, C = _generate_body()
    parent_interframe = ReferenceFrame('int_frame')
    parent_interframe.orient_axis(N, N.z, 1)
    joint = PinJoint('J', P, C, parent_interframe=parent_interframe)
    assert joint.parent_interframe == parent_interframe
    assert joint.parent_interframe.ang_vel_in(N) == 0
    assert joint.child_interframe == A
    # Test invalid type
    N, A, P, C = _generate_body()
    raises(TypeError, lambda: PinJoint('J', P, C, parent_interframe=N.x))
    # Test time varying orientations
    q = dynamicsymbols('q')
    N, A, P, C = _generate_body()
    parent_interframe = ReferenceFrame('int_frame')
    parent_interframe.orient_axis(N, N.z, q)
    raises(ValueError,
           lambda: PinJoint('J', P, C, parent_interframe=parent_interframe))
    # Test undefined frame
    N, A, P, C = _generate_body()
    child_interframe = ReferenceFrame('int_frame')
    child_interframe.orient_axis(N, N.z, 1)  # Defined with respect to parent
    raises(ValueError,
           lambda: PinJoint('J', P, C, child_interframe=child_interframe))


def test_slidingjoint():
    _, _, P, C = _generate_body()
    x, v = dynamicsymbols('x_S, v_S')
    S = PrismaticJoint('S', P, C)
    assert S.name == 'S'
    assert S.parent == P
    assert S.child == C
    assert S.coordinates == [x]
    assert S.speeds == [v]
    assert S.kdes == [v - x.diff(t)]
    assert S.joint_axis == P.frame.x
    assert S.child_point.pos_from(C.masscenter) == Vector(0)
    assert S.parent_point.pos_from(P.masscenter) == Vector(0)
    assert S.parent_point.pos_from(S.child_point) == - x * P.frame.x
    assert P.masscenter.pos_from(C.masscenter) == - x * P.frame.x
    assert C.masscenter.vel(P.frame) == v * P.frame.x
    assert P.ang_vel_in(C) == 0
    assert C.ang_vel_in(P) == 0
    assert S.__str__() == 'PrismaticJoint: S  parent: P  child: C'

    N, A, P, C = _generate_body()
    l, m = symbols('l m')
    Pint = ReferenceFrame('P_int')
    Pint.orient_axis(P.frame, P.y, pi / 2)
    S = PrismaticJoint('S', P, C, parent_point=l * P.frame.x,
                       child_point=m * C.frame.y, joint_axis=P.frame.z,
                       parent_interframe=Pint)

    assert S.joint_axis == P.frame.z
    assert S.child_point.pos_from(C.masscenter) == m * C.frame.y
    assert S.parent_point.pos_from(P.masscenter) == l * P.frame.x
    assert S.parent_point.pos_from(S.child_point) == - x * P.frame.z
    assert P.masscenter.pos_from(C.masscenter) == - l*N.x - x*N.z + m*A.y
    assert C.masscenter.vel(P.frame) == v * P.frame.z
    assert C.ang_vel_in(P) == 0
    assert P.ang_vel_in(C) == 0

    _, _, P, C = _generate_body()
    Pint = ReferenceFrame('P_int')
    Pint.orient_axis(P.frame, P.y, pi / 2)
    S = PrismaticJoint('S', P, C, parent_point=l * P.frame.z,
                       child_point=m * C.frame.x, joint_axis=P.frame.z,
                       parent_interframe=Pint)
    assert S.joint_axis == P.frame.z
    assert S.child_point.pos_from(C.masscenter) == m * C.frame.x
    assert S.parent_point.pos_from(P.masscenter) == l * P.frame.z
    assert S.parent_point.pos_from(S.child_point) == - x * P.frame.z
    assert P.masscenter.pos_from(C.masscenter) == (-l - x)*P.frame.z + m*C.frame.x
    assert C.masscenter.vel(P.frame) == v * P.frame.z
    assert C.ang_vel_in(P) == 0
    assert P.ang_vel_in(C) == 0


def test_slidingjoint_arbitrary_axis():
    x, v = dynamicsymbols('x_S, v_S')

    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PrismaticJoint('S', P, C, child_axis=-A.x)

    assert (-A.x).angle_between(N.x) == 0
    assert -A.x.express(N) == N.x
    assert A.dcm(N) == Matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    assert C.masscenter.pos_from(P.masscenter) == x * N.x
    assert C.masscenter.pos_from(P.masscenter).express(A).simplify() == -x * A.x
    assert C.masscenter.vel(N) == v * N.x
    assert C.masscenter.vel(N).express(A) == -v * A.x
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    #When axes are different and parent joint is at masscenter but child joint is at a unit vector from
    #child masscenter.
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PrismaticJoint('S', P, C, child_axis=A.y, child_point=A.x)

    assert A.y.angle_between(N.x) == 0 #Axis are aligned
    assert A.y.express(N) == N.x
    assert A.dcm(N) == Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert C.masscenter.vel(N) == v * N.x
    assert C.masscenter.vel(N).express(A) == v * A.y
    assert C.masscenter.pos_from(P.masscenter) == x*N.x - A.x
    assert C.masscenter.pos_from(P.masscenter).express(N).simplify() == x*N.x + N.y
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    #Similar to previous case but wrt parent body
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PrismaticJoint('S', P, C, parent_axis=N.y, parent_point=N.x)

    assert N.y.angle_between(A.x) == 0 #Axis are aligned
    assert N.y.express(A) ==  A.x
    assert A.dcm(N) == Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    assert C.masscenter.vel(N) == v * N.y
    assert C.masscenter.vel(N).express(A) == v * A.x
    assert C.masscenter.pos_from(P.masscenter) == N.x + x*N.y
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    #Both joint pos is defined but different axes
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PrismaticJoint('S', P, C, parent_point=N.x, child_point=A.x,
                       child_axis=A.x+A.y)
    assert N.x.angle_between(A.x + A.y) == 0 #Axis are aligned
    assert (A.x + A.y).express(N) == sqrt(2)*N.x
    assert A.dcm(N) == Matrix([[sqrt(2)/2, -sqrt(2)/2, 0], [sqrt(2)/2, sqrt(2)/2, 0], [0, 0, 1]])
    assert C.masscenter.pos_from(P.masscenter) == (x + 1)*N.x - A.x
    assert C.masscenter.pos_from(P.masscenter).express(N) == (x - sqrt(2)/2 + 1)*N.x + sqrt(2)/2*N.y
    assert C.masscenter.vel(N).express(A) == v * (A.x + A.y)/sqrt(2)
    assert C.masscenter.vel(N) == v*N.x
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        PrismaticJoint('S', P, C, parent_point=N.x, child_point=A.x,
                       child_axis=A.x+A.y-A.z)
    assert N.x.angle_between(A.x + A.y - A.z) == 0 #Axis are aligned
    assert (A.x + A.y - A.z).express(N) == sqrt(3)*N.x
    assert _simplify_matrix(A.dcm(N)) == Matrix([[sqrt(3)/3, -sqrt(3)/3, sqrt(3)/3],
                                                 [sqrt(3)/3, sqrt(3)/6 + S(1)/2, S(1)/2 - sqrt(3)/6],
                                                 [-sqrt(3)/3, S(1)/2 - sqrt(3)/6, sqrt(3)/6 + S(1)/2]])
    assert C.masscenter.pos_from(P.masscenter) == (x + 1)*N.x - A.x
    assert C.masscenter.pos_from(P.masscenter).express(N) == \
        (x - sqrt(3)/3 + 1)*N.x + sqrt(3)/3*N.y - sqrt(3)/3*N.z
    assert C.masscenter.vel(N) == v*N.x
    assert C.masscenter.vel(N).express(A) == sqrt(3)*v/3*A.x + sqrt(3)*v/3*A.y - sqrt(3)*v/3*A.z
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0

    N, A, P, C = _generate_body()
    m, n = symbols('m n')
    with warns_deprecated_sympy():
        PrismaticJoint('S', P, C, parent_point=m*N.x, child_point=n*A.x,
                       child_axis=A.x+A.y-A.z, parent_axis=N.x-N.y+N.z)
    assert (N.x-N.y+N.z).angle_between(A.x+A.y-A.z) == 0 #Axis are aligned
    assert (A.x+A.y-A.z).express(N) == N.x - N.y + N.z
    assert _simplify_matrix(A.dcm(N)) == Matrix([[-S(1)/3, -S(2)/3, S(2)/3],
                                                 [S(2)/3, S(1)/3, S(2)/3],
                                                 [-S(2)/3, S(2)/3, S(1)/3]])
    assert C.masscenter.pos_from(P.masscenter) == \
        (m + sqrt(3)*x/3)*N.x - sqrt(3)*x/3*N.y + sqrt(3)*x/3*N.z - n*A.x
    assert C.masscenter.pos_from(P.masscenter).express(N) == \
        (m + n/3 + sqrt(3)*x/3)*N.x + (2*n/3 - sqrt(3)*x/3)*N.y + (-2*n/3 + sqrt(3)*x/3)*N.z
    assert C.masscenter.vel(N) == sqrt(3)*v/3*N.x - sqrt(3)*v/3*N.y + sqrt(3)*v/3*N.z
    assert C.masscenter.vel(N).express(A) == sqrt(3)*v/3*A.x + sqrt(3)*v/3*A.y - sqrt(3)*v/3*A.z
    assert A.ang_vel_in(N) == 0
    assert N.ang_vel_in(A) == 0


def test_deprecated_joint_pos():
    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        pin = PinJoint('J', P, C, parent_joint_pos=N.x + N.y,
                       child_joint_pos=C.y - C.z)
    assert pin.parent_point.pos_from(P.masscenter) == N.x + N.y
    assert pin.child_point.pos_from(C.masscenter) == C.y - C.z

    N, A, P, C = _generate_body()
    with warns_deprecated_sympy():
        slider = PrismaticJoint('J', P, C, parent_joint_pos=N.z + N.y,
                                child_joint_pos=C.y - C.x)
    assert slider.parent_point.pos_from(P.masscenter) == N.z + N.y
    assert slider.child_point.pos_from(C.masscenter) == C.y - C.x
