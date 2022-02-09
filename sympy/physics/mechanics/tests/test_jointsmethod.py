from sympy.core.function import expand
from sympy.core.symbol import symbols
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.matrices.dense import Matrix
from sympy.simplify.trigsimp import trigsimp
from sympy.physics.mechanics import (PinJoint, JointsMethod, Body, KanesMethod,
                                    PrismaticJoint, LagrangesMethod, inertia)
from sympy.physics.vector import dynamicsymbols, ReferenceFrame
from sympy.testing.pytest import raises


t = dynamicsymbols._t # type: ignore


def test_jointsmethod():
    P = Body('P')
    C = Body('C')
    Pin = PinJoint('P1', P, C)
    C_ixx, g = symbols('C_ixx g')
    theta, omega = dynamicsymbols('theta_P1, omega_P1')
    P.apply_force(g*P.y)
    method = JointsMethod(P, Pin)
    assert method.frame == P.frame
    assert method.bodies == [C, P]
    assert method.loads == [(P.masscenter, g*P.frame.y)]
    assert method.q == [theta]
    assert method.u == [omega]
    assert method.kdes == [omega - theta.diff()]
    soln = method.form_eoms()
    assert soln == Matrix([[-C_ixx*omega.diff()]])
    assert method.forcing_full == Matrix([[omega], [0]])
    assert method.mass_matrix_full == Matrix([[1, 0], [0, C_ixx]])
    assert isinstance(method.method, KanesMethod)

def test_jointmethod_duplicate_coordinates_speeds():
    P = Body('P')
    C = Body('C')
    T = Body('T')
    q, u = dynamicsymbols('q u')
    P1 = PinJoint('P1', P, C, q)
    P2 = PrismaticJoint('P2', C, T, q)
    raises(ValueError, lambda: JointsMethod(P, P1, P2))

    P1 = PinJoint('P1', P, C, speeds=u)
    P2 = PrismaticJoint('P2', C, T, speeds=u)
    raises(ValueError, lambda: JointsMethod(P, P1, P2))

    P1 = PinJoint('P1', P, C, q, u)
    P2 = PrismaticJoint('P2', C, T, q, u)
    raises(ValueError, lambda: JointsMethod(P, P1, P2))

def test_complete_simple_double_pendulum():
    q1, q2 = dynamicsymbols('q1 q2')
    u1, u2 = dynamicsymbols('u1 u2')
    m, l, g = symbols('m l g')
    C = Body('C')  # ceiling
    PartP = Body('P', mass=m)
    PartR = Body('R', mass=m)

    J1 = PinJoint('J1', C, PartP, speeds=u1, coordinates=q1,
                  child_joint_pos=-l*PartP.x, parent_axis=C.z,
                  child_axis=PartP.z)
    J2 = PinJoint('J2', PartP, PartR, speeds=u2, coordinates=q2,
                  child_joint_pos=-l*PartR.x, parent_axis=PartP.z,
                  child_axis=PartR.z)

    PartP.apply_force(m*g*C.x)
    PartR.apply_force(m*g*C.x)

    method = JointsMethod(C, J1, J2)
    method.form_eoms()

    assert expand(method.mass_matrix_full) == Matrix([[1, 0, 0, 0],
                                                      [0, 1, 0, 0],
                                                      [0, 0, 2*l**2*m*cos(q2) + 3*l**2*m, l**2*m*cos(q2) + l**2*m],
                                                      [0, 0, l**2*m*cos(q2) + l**2*m, l**2*m]])
    assert trigsimp(method.forcing_full) == trigsimp(Matrix([[u1], [u2], [-g*l*m*(sin(q1 + q2) + sin(q1)) -
                                           g*l*m*sin(q1) + l**2*m*(2*u1 + u2)*u2*sin(q2)],
                                          [-g*l*m*sin(q1 + q2) - l**2*m*u1**2*sin(q2)]]))

def test_two_dof_joints():
    q1, q2, u1, u2 = dynamicsymbols('q1 q2 u1 u2')
    m, c1, c2, k1, k2 = symbols('m c1 c2 k1 k2')
    W = Body('W')
    B1 = Body('B1', mass=m)
    B2 = Body('B2', mass=m)
    J1 = PrismaticJoint('J1', W, B1, coordinates=q1, speeds=u1)
    J2 = PrismaticJoint('J2', B1, B2, coordinates=q2, speeds=u2)
    W.apply_force(k1*q1*W.x, reaction_body=B1)
    W.apply_force(c1*u1*W.x, reaction_body=B1)
    B1.apply_force(k2*q2*W.x, reaction_body=B2)
    B1.apply_force(c2*u2*W.x, reaction_body=B2)
    method = JointsMethod(W, J1, J2)
    method.form_eoms()
    MM = method.mass_matrix
    forcing = method.forcing
    rhs = MM.LUsolve(forcing)
    assert expand(rhs[0]) == expand((-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2)/m)
    assert expand(rhs[1]) == expand((k1 * q1 + c1 * u1 - 2 * k2 * q2 - 2 *
                                    c2 * u2) / m)

def test_simple_pedulum():
    l, m, g = symbols('l m g')
    C = Body('C')
    b = Body('b', mass=m)
    q = dynamicsymbols('q')
    P = PinJoint('P', C, b, speeds=q.diff(t), coordinates=q, child_joint_pos = -l*b.x,
                    parent_axis=C.z, child_axis=b.z)
    b.potential_energy = - m * g * l * cos(q)
    method = JointsMethod(C, P)
    method.form_eoms(LagrangesMethod)
    rhs = method.rhs()
    assert rhs[1] == -g*sin(q)/l

def test_chaos_pendulum():
    #https://www.pydy.org/examples/chaos_pendulum.html
    mA, mB, lA, lB, IAxx, IBxx, IByy, IBzz, g = symbols('mA, mB, lA, lB, IAxx, IBxx, IByy, IBzz, g')
    theta, phi, omega, alpha = dynamicsymbols('theta phi omega alpha')

    A = ReferenceFrame('A')
    B = ReferenceFrame('B')

    rod = Body('rod', mass=mA, frame=A, central_inertia=inertia(A, IAxx, IAxx, 0))
    plate = Body('plate', mass=mB, frame=B, central_inertia=inertia(B, IBxx, IByy, IBzz))
    C = Body('C')
    J1 = PinJoint('J1', C, rod, coordinates=theta, speeds=omega,
                  child_joint_pos=-lA*rod.z, parent_axis=C.y, child_axis=rod.y)
    J2 = PinJoint('J2', rod, plate, coordinates=phi, speeds=alpha,
                  parent_joint_pos=(lB-lA)*rod.z, parent_axis=rod.z, child_axis=plate.z)

    rod.apply_force(mA*g*C.z)
    plate.apply_force(mB*g*C.z)

    method = JointsMethod(C, J1, J2)
    method.form_eoms()

    MM = method.mass_matrix
    forcing = method.forcing
    rhs = MM.LUsolve(forcing)
    xd = (-2 * IBxx * alpha * omega * sin(phi) * cos(phi) + 2 * IByy * alpha * omega * sin(phi) *
            cos(phi) - g * lA * mA * sin(theta) - g * lB * mB * sin(theta)) / (IAxx + IBxx *
                sin(phi)**2 + IByy * cos(phi)**2 + lA**2 * mA + lB**2 * mB)
    assert (rhs[0] - xd).simplify() == 0
    xd = (IBxx - IByy) * omega**2 * sin(phi) * cos(phi) / IBzz
    assert (rhs[1] - xd).simplify() == 0
