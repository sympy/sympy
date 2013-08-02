from sympy import symbols
from sympy.physics.mechanics import Particle, MovingRefFrame, dynamicsymbols


def test_particle():
    """ Tests for Particle class """

    R0 = MovingRefFrame('R0')
    a, b, c = symbols('a b c')
    v1, v2, v3 = symbols('v1 v2 v3')
    v = a*R0.x + b*R0.y + c*R0.z
    V = v1*R0.x + v2*R0.y + v3*R0.z
    m, p = symbols('m p')
    P0 = Particle('P0', R0, m, p)
    assert P0.frame = R0
    assert P0.mass == m
    assert P0.potential_energy == p
    assert P0.linear_momentum(R0) == 0
    assert P0.angular_momentum(V, R0) == 0
    R1 = MovingRefFrame('R1', trans_vel = v, parentframe = R0)
    assert P0.linear_momentum(R1) == -m*v
    P1 = Particle('P1', R1, m, p)
    assert P1.pos_vector_wrt(P0) == P1.pos_vector_wrt(R0)
    assert P1.trans_vel_wrt(P0) == P1.trans_vel_wrt(R0)
    assert P1.trans_acc_wrt(P0) == P1.trans_acc_wrt(R0)
    assert P1.linear_momentum(R0) == m*v
    assert P1.angular_momentum(3 * R0.x, R0) == \
           (-b*m*(c*t - v3) + c*m*(b*t - v2))*R0.x + \
           (a*m*(c*t - v3) - c*m*(a*t - v1))*R0.y + \
           (-a*m*(b*t - v2) + b*m*(a*t - v1))*R0.z
    assert P1.kinetic_energy(R1) == 0
    assert P1.total_energy(R1) == p
    assert P1.kinetic_energy(R0) == a**2*m/2 + b**2*m/2 + c**2*m/2
    assert P1.total_energy(R0) == p + a**2*m/2 + b**2*m/2 + c**2*m/2
    
