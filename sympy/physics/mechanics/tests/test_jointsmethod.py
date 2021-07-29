from sympy import symbols, Matrix, cos, sin, expand
from sympy.physics.mechanics import (PinJoint, JointsMethod, Body, KanesMethod,
                                    PrismaticJoint, LagrangesMethod)
from sympy.physics.vector import dynamicsymbols, ReferenceFrame


t = dynamicsymbols._t


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

def test_complete_simple_double_pendulum():
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

    method = JointsMethod(C, J1, J2)
    method.form_eoms()
    assert expand(method.mass_matrix_full) == Matrix([[1, 0, 0, 0], [0, 1, 0, 0],
                                        [0, 0, 2*l**2*m*cos(q2) + 3*l**2*m, l**2*m*cos(q2) + l**2*m],
                                        [0, 0, l**2*m*cos(q2) + l**2*m, l**2*m]])
    assert expand(method.forcing_full) == Matrix([[u1], [u2], [2*l**2*m*u1*u2*sin(q2) +
                                    l**2*m*u2**2*sin(q2)], [-l**2*m*u1**2*sin(q2)]])


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
