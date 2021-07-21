from sympy import symbols, Matrix, cos, sin
from sympy.physics.mechanics import PinJoint, JointsMethod, Body, KanesMethod
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
    (fr, frstar) = method.form_eoms()
    mm = method.mass_matrix_full
    fo = method.forcing_full
    qudots = mm.inv() * fo
    qudots.simplify()
    qudots == Matrix([[o1], [o2], [(l2*(-(o1 + o2)**2 + o1**2) - (l1*cos(t2) + l2)*o1**2)*sin(t2)/(l1*(cos(t2)**2 - 2))], [(l2*(l1*cos(t2) + l2)*((o1 + o2)**2 - o1**2) + (2*l1**2 + 2*l1*l2*cos(t2) + l2**2)*o1**2)*sin(t2)/(l1*l2*(cos(t2)**2 - 2))]])
