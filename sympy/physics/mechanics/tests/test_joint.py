from sympy.core.symbol import symbols
from sympy import Symbol, sin, cos, Matrix
from sympy.physics.vector import Vector, cross, ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols, Body, PinJoint, inertia
from sympy.physics.mechanics.joint import Joint
from sympy.testing.pytest import raises

t = dynamicsymbols._t

def test_Joint():
    parent = Body('parent')
    child = Body('child')
    raises(TypeError, lambda: Joint('J', parent, child))


def test_pinjoint():
    P = Body('P')
    C = Body('C')
    Pj = PinJoint('J', P, C)
    theta, omega = dynamicsymbols('J_theta, J_omega')

    assert Pj._name == 'J'
    assert Pj.parent() == P
    assert Pj.child() == C
    assert Pj.coordinates() == [theta]
    assert Pj.speeds() == [omega]
    assert Pj.kdes() == [omega - theta.diff(t)]
    assert Pj._parent_axis == P.frame.x
    assert Pj._child_joint.pos_from(C.masscenter) == Vector(0)
    assert Pj._parent_joint.pos_from(P.masscenter) == Vector(0)
    assert Pj._parent_joint.pos_from(Pj._child_joint) == Vector(0)

    l, m = symbols('l m')
    J1 = PinJoint('J1', P, C, parent_joint_pos= l * P.frame.x, child_joint_pos= m * C.frame.y,
        parent_axis = P.frame.z)

    assert J1._parent_axis == P.frame.z
    assert J1._child_joint.pos_from(C.masscenter) == m * C.frame.y
    assert J1._parent_joint.pos_from(P.masscenter) == l * P.frame.x
    assert J1._parent_joint.pos_from(J1._child_joint) == Vector(0)

def test_pin_joint_double_pendulum():
    q1, q2 = dynamicsymbols('q1 q2')
    u1, u2 = dynamicsymbols('u1 u2')
    m = symbols('m')
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    B = ReferenceFrame('B')
    C = Body('C', frame=N) #ceiling
    PartP = Body('P', frame=A, mass=m)
    PartR = Body('R', frame=B, mass=m)

    J1 = PinJoint('J1', C, PartP, speeds=u1, coordinates=q1, parent_axis=C.frame.z)
    J2 = PinJoint('J2', PartP, PartR, speeds=u2, coordinates=q2, parent_axis=PartP.frame.z)

    #Check orientation
    assert N.dcm(A) == Matrix([[cos(q1), -sin(q1), 0], [sin(q1), cos(q1), 0], [0, 0, 1]])
    assert A.dcm(B) == Matrix([[cos(q2), -sin(q2), 0], [sin(q2), cos(q2), 0], [0, 0, 1]])
    assert N.dcm(B).simplify() == Matrix([[cos(q1 + q2), -sin(q1 + q2), 0], \
        [sin(q1 + q2), cos(q1 + q2), 0], [0, 0, 1]])

    #Check Angular Velocity
    assert A.ang_vel_in(N) == u1 * N.z
    assert B.ang_vel_in(A) == u2 * A.z
    assert B.ang_vel_in(N) == u1 * N.z + u2 * A.z

    #Check kde
    assert J1.kdes == [u1 - q1.diff(t)]
    assert J1.kdes() == [u2 - q2.diff(t)]

    #Check Linear Velocity
    assert J1._parent_joint.vel(N) == 0
    # Problem here

def test_pin_joint_chaos_pendulum():
    mA, mB, lA, lB = symbols('mA, mB, lA, lB')
    theta, phi, omega, alpha = dynamicsymbols('theta phi omega alpha')
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    B = ReferenceFrame('B')
    No = Point('No')
    Ao = Point('Ao')
    Bo = Point('Bo')
    Ao.set_pos(No, lA * A.z)
    Bo.set_pos(No, lB * A.z)
    rod = Body('rod', Ao, mA, A)
    plate = Body('plate', Bo, mB, B)
    C = Body('C', No, frame=N)
    J1 = PinJoint('J1', C, rod, coordinates=theta, speeds=omega, parent_axis=N.y)
    J2 = PinJoint('J2',rod, plate, coordinates=phi, speeds=alpha, parent_axis=A.z)

    #Check orientation
    assert A.dcm(N) == Matrix([[cos(theta), 0, -sin(theta)], [0, 1, 0], [sin(theta), 0, cos(theta)]])
    assert A.dcm(B) == Matrix([[cos(phi), -sin(phi), 0], [sin(phi), cos(phi), 0], [0, 0, 1]])
    assert B.dcm(N) == Matrix([[cos(phi)*cos(theta), sin(phi), -sin(theta)*cos(phi)], \
        [-sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta)], [sin(theta), 0, cos(theta)]])

    #Check Angular Velocity
    assert A.ang_vel_in(N) == omega * N.y
    assert A.ang_vel_in(B) == - alpha * A.z
    assert N.ang_vel_in(B) == - omega * N.y - alpha * A.z

    #Check kde
    assert J1.kdes() == [omega - theta.diff(t)]
    assert J2.kdes() == [alpha - phi.diff(t)]

    #Check Linear Velocities
    assert J1._parent_joint.vel(N) == 0
    # Incorrect result in linear velocity
    # J1._child_joint.vel(A)
