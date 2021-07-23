from sympy import symbols, Matrix, cos, sin, expand
from sympy.physics.mechanics import PinJoint, JointsMethod, Body, KanesMethod, PrismaticJoint
from sympy.physics.vector import dynamicsymbols

def test_jointsmethod():
    P = Body('P')
    C = Body('C')
    Pin = PinJoint('P1', P, C)
    C_ixx, g = symbols('C_ixx g')
    theta, omega = dynamicsymbols('theta_P1, omega_P1')
    P.apply_force(g*P.y)
    method = JointsMethod(P, Pin)
    assert method.frame == P.frame
    assert method.bodylist == [C, P]
    assert method.loadlist == [(P.masscenter, g*P.frame.y)]
    assert method.q == [theta]
    assert method.u == [omega]
    assert method.kdes == [omega - theta.diff()]
    fr, frstar = method.form_eoms()
    assert fr == Matrix([[0]])
    assert frstar == Matrix([[-C_ixx*omega.diff()]])
    assert method.forcing_full == Matrix([[omega], [0]])
    assert method.mass_matrix_full == Matrix([[1, 0], [0, C_ixx]])
    assert isinstance(method.method, KanesMethod)

def test_complete_double_pendulum():
    l1, l2, m = symbols('l1 l2 m')
    o1, o2 = dynamicsymbols('omega_P1, omega_P2')
    t1, t2 = dynamicsymbols('theta_P1, theta_P2')
    ceiling = Body('C')
    upper_bob = Body('U', mass=m)
    lower_bob = Body('L', mass=m)
    ceiling_joint = PinJoint('P1', ceiling, upper_bob,child_joint_pos=-l1*upper_bob.frame.x,
        parent_axis=ceiling.frame.z, child_axis=upper_bob.frame.z)
    pendulum_joint = PinJoint('P2', upper_bob, lower_bob,child_joint_pos=-l2*lower_bob.frame.x,
        parent_axis=upper_bob.frame.z,child_axis=lower_bob.frame.z)
    method = JointsMethod(ceiling, ceiling_joint, pendulum_joint)
    method.form_eoms()
    mm = method.mass_matrix_full
    fo = method.forcing_full
    qudots = mm.inv() * fo
    qudots.simplify()
    qudots == Matrix([[o1], [o2], [(l2*(-(o1 + o2)**2 + o1**2) - (l1*cos(t2) + l2)*o1**2)*sin(t2)/(l1*(cos(t2)**2 - 2))], [(l2*(l1*cos(t2) + l2)*((o1 + o2)**2 - o1**2) + (2*l1**2 + 2*l1*l2*cos(t2) + l2**2)*o1**2)*sin(t2)/(l1*l2*(cos(t2)**2 - 2))]])

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
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand((-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2)/m)
    assert expand(rhs[1]) == expand((k1 * q1 + c1 * u1 - 2 * k2 * q2 - 2 *
                                    c2 * u2) / m)
