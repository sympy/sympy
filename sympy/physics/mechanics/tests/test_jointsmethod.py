from sympy import symbols, Matrix
from sympy.physics.mechanics import PinJoint, JointsMethod, Body
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
    fr, frstar = method.kanes_equations()
    assert fr == Matrix([[0]])
    assert frstar == Matrix([[-C_ixx*omega.diff()]])
    assert method.forcing_full == Matrix([[omega], [0]])
    assert method.mass_matrix_full == Matrix([[1, 0], [0, C_ixx]])
    assert method.fr == fr
    assert method.frstar == frstar
