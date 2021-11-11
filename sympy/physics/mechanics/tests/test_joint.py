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
from sympy.physics.vector import Vector, ReferenceFrame
from sympy.testing.pytest import raises

t = dynamicsymbols._t # type: ignore


def _generate_body():
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    P = Body('P', frame=N)
    C = Body('C', frame=A)
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
    assert Pj.parent_axis == P.frame.x
    assert Pj.child_point.pos_from(C.masscenter) == Vector(0)
    assert Pj.parent_point.pos_from(P.masscenter) == Vector(0)
    assert Pj.parent_point.pos_from(Pj._child_point) == Vector(0)
    assert C.masscenter.pos_from(P.masscenter) == Vector(0)
    assert Pj.__str__() == 'PinJoint: J  parent: P  child: C'

    P1 = Body('P1')
    C1 = Body('C1')
    J1 = PinJoint('J1', P1, C1, parent_joint_pos=l*P1.frame.x,
                  child_joint_pos=m*C1.frame.y, parent_axis=P1.frame.z)

    assert J1._parent_axis == P1.frame.z
    assert J1._child_point.pos_from(C1.masscenter) == m * C1.frame.y
    assert J1._parent_point.pos_from(P1.masscenter) == l * P1.frame.x
    assert J1._parent_point.pos_from(J1._child_point) == Vector(0)
    assert (P1.masscenter.pos_from(C1.masscenter) ==
            -l*P1.frame.x + m*C1.frame.y)


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
                  child_joint_pos=-l*A.x, parent_axis=C.frame.z,
                  child_axis=PartP.frame.z)
    J2 = PinJoint('J2', PartP, PartR, speeds=u2, coordinates=q2,
                  child_joint_pos=-l*B.x, parent_axis=PartP.frame.z,
                  child_axis=PartR.frame.z)

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
                  child_joint_pos=lA*A.z, parent_axis=N.y, child_axis=A.y)
    J2 = PinJoint('J2', rod, plate, coordinates=phi, speeds=alpha,
                  parent_joint_pos=lC*A.z, parent_axis=A.z, child_axis=B.z)

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


def test_pinjoint_arbitrary_axis():
    theta, omega = dynamicsymbols('theta_J, omega_J')

    # When the bodies are attached though masscenters but axess are opposite.
    N, A, P, C = _generate_body()
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
    PinJoint('J', P, C, child_axis=A.y, child_joint_pos=A.x)

    assert A.y.angle_between(N.x) == 0  # Axis are aligned
    assert A.y.express(N) == N.x
    assert A.dcm(N) == Matrix([[0, -cos(theta), -sin(theta)],
                               [1, 0, 0],
                               [0, -sin(theta), cos(theta)]])
    assert A.ang_vel_in(N) == omega*N.x
    assert A.ang_vel_in(N).express(A) == omega * A.y
    assert A.ang_vel_in(N).magnitude() == sqrt(omega**2)
    angle = A.ang_vel_in(N).angle_between(A.y)
    assert angle.xreplace({omega: 1}) == 0
    assert C.masscenter.vel(N) == omega*A.z
    assert C.masscenter.pos_from(P.masscenter) == -A.x
    assert (C.masscenter.pos_from(P.masscenter).express(N).simplify() ==
            cos(theta)*N.y + sin(theta)*N.z)
    assert C.masscenter.vel(N).angle_between(A.x) == pi/2

    # Similar to previous case but wrt parent body
    N, A, P, C = _generate_body()
    PinJoint('J', P, C, parent_axis=N.y, parent_joint_pos=N.x)

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
    assert C.masscenter.vel(N).simplify() == - omega*N.z
    assert C.masscenter.pos_from(P.masscenter) == N.x

    # Both joint pos id defined but different axes
    N, A, P, C = _generate_body()
    PinJoint('J', P, C, parent_joint_pos=N.x, child_joint_pos=A.x,
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
    PinJoint('J', P, C, parent_joint_pos=N.x, child_joint_pos=A.x,
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
    PinJoint('J', P, C, parent_joint_pos=m*N.x, child_joint_pos=n*A.x,
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
            (m*omega*N.y + m*omega*N.z + n*omega*A.y + n*omega*A.z)/sqrt(3))
    assert C.masscenter.pos_from(P.masscenter) == m*N.x - n*A.x
    assert (C.masscenter.pos_from(P.masscenter).express(N).simplify() ==
            (m + n*(2*cos(theta) - 1)/3)*N.x + n*(2*sin(theta + pi/6) +
            1)/3*N.y - n*(2*cos(theta + pi/3) + 1)/3*N.z)
    assert (C.masscenter.vel(N).express(N).simplify() ==
            -2*n*omega*sin(theta)/3*N.x + (sqrt(3)*m +
                                           2*n*cos(theta + pi/6))*omega/3*N.y
            + (sqrt(3)*m + 2*n*sin(theta + pi/3))*omega/3*N.z)
    assert expand_mul(C.masscenter.vel(N).angle_between(m*N.x - n*A.x)) == pi/2


def test_pinjoint_pi():
    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, child_axis=-C.frame.x)
    assert J._generate_vector() == P.frame.z

    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, parent_axis=P.frame.y, child_axis=-C.frame.y)
    assert J._generate_vector() == P.frame.x

    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, parent_axis=P.frame.z, child_axis=-C.frame.z)
    assert J._generate_vector() == P.frame.y

    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, parent_axis=P.frame.x+P.frame.y, child_axis=-C.frame.y-C.frame.x)
    assert J._generate_vector() == P.frame.z

    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, parent_axis=P.frame.y+P.frame.z, child_axis=-C.frame.y-C.frame.z)
    assert J._generate_vector() == P.frame.x

    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, parent_axis=P.frame.x+P.frame.z, child_axis=-C.frame.z-C.frame.x)
    assert J._generate_vector() == P.frame.y

    _, _, P, C = _generate_body()
    J = PinJoint('J', P, C, parent_axis=P.frame.x+P.frame.y+P.frame.z,
                child_axis=-C.frame.x-C.frame.y-C.frame.z)
    assert J._generate_vector() == P.frame.y - P.frame.z


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
    assert S.parent_axis == P.frame.x
    assert S.child_axis == C.frame.x
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
    S = PrismaticJoint('S', P, C, parent_joint_pos= l * P.frame.x,
                    child_joint_pos= m * C.frame.y, parent_axis = P.frame.z)

    assert S.parent_axis == P.frame.z
    assert S.child_point.pos_from(C.masscenter) == m * C.frame.y
    assert S.parent_point.pos_from(P.masscenter) == l * P.frame.x
    assert S.parent_point.pos_from(S.child_point) == - x * P.frame.z
    assert P.masscenter.pos_from(C.masscenter) == - l*N.x - x*N.z + m*A.y
    assert C.masscenter.vel(P.frame) == v * P.frame.z
    assert C.ang_vel_in(P) == 0
    assert P.ang_vel_in(C) == 0

    _, _, P, C = _generate_body()
    S = PrismaticJoint('S', P, C, parent_joint_pos= l * P.frame.z,
                    child_joint_pos= m * C.frame.x, parent_axis = P.frame.z)
    assert S.parent_axis == P.frame.z
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
    PrismaticJoint('S', P, C, child_axis=A.y, child_joint_pos=A.x)

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
    PrismaticJoint('S', P, C, parent_axis=N.y, parent_joint_pos=N.x)

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
    PrismaticJoint('S', P, C, parent_joint_pos=N.x, child_joint_pos=A.x, child_axis=A.x+A.y)
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
    PrismaticJoint('S', P, C, parent_joint_pos=N.x, child_joint_pos=A.x, child_axis=A.x+A.y-A.z)
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
    PrismaticJoint('S', P, C, parent_joint_pos=m*N.x, child_joint_pos=n*A.x,
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
