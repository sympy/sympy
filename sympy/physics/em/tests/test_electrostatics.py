from sympy import symbols
from sympy.physics.em import ParticleCharge, gradient, scalar_potential, \
     charge_density, charge_assembly_energy
from sympy.physics.mechanics import ReferenceFrame, get_motion_acc, \
     dynamicsymbols
from sympy.utilities.pytest import raises

R = ReferenceFrame('R')


def test_particle_charge():
    m, q = symbols('m q')
    P = ParticleCharge('P', m, q)
    assert raises(ValueError, P.electrostatic_potential(R))
    assert P.charge == q
    P.set_motion(R, pos_vector = 0)
    assert gradient(P.electrostatic_potential(R), R) == \
           -P.electrostatic_field(R)
    assert P.electrostatic_potential(R) == \
           k*q*(R[0]**2 + R[1]**2 + R[2]**2)**(-0.5)
    assert P.electrostatic_potential(R, R.x+R.y+R.z) == 3**(-0.5)*k*q
    assert P.electrostatic_field(R) == \
           - 1.0*R[0]*k*q*(R[0]**2 + R[1]**2 + R[2]**2)**(-1.5)*R.x - \
           1.0*R[1]*k*q*(R[0]**2 + R[1]**2 + R[2]**2)**(-1.5)*R.y - \
           1.0*R[2]*k*q*(R[0]**2 + R[1]**2 + R[2]**2)**(-1.5)*R.z
    assert P.electrostatic_field(R, R.x+R.y+R.z) == - 1.0*3**(-1.5)*k*q*R.x - \
           1.0*3**(-1.5)*k*q*R.y - 1.0*3**(-1.5)*k*q*R.z
    assert P.electrostatic_force(R.x) == q*R.x
    assert P.electrostatic_force(R.x, R.y, R.z) == \
           P.electrostatic_force(R.x+R.y+R.z) == q*(R.x+R.y+R.z)
    P.set_motion(R, pos_vector=R.x)
    field1 = P.electrostatic_field(R, 0)
    pot1 = P.electrostatic_potential(R, 0)
    P.set_motion(R, pos_vector = -R.x)
    field2 = P.electrostatic_field(R, 0)
    pot2 = P.electrostatic_potential(R, 0)
    assert pot1 == pot2
    assert field1 + field2 == 0


def test_field_and_potential():
    potential = 2*R[0]**2*R[1]*R[2]
    field = -gradient(scalar_field, R)
    assert electrostatic_potential(field, frame) == potential
    assert electrostatic_field(potential, R) == field


def test_charge_density():
    a,b,c = symbols('a b c')
    field = a*R[0]*R.x + b*R[1]*R.y + c*R[2]*R.z
    assert charge_density(field, R) == (a + b + c)/(4*k*pi)
    potential1 = R[0]*R[1]*R[2]
    assert charge_density(potential1, R) == 0
    potential2 = potential1**2
    assert charge_density(potential2, R) == (-2*R[0]**2*R[1]**2 - \
                                             2*R[0]**2*R[2]**2 - \
                                             2*R[1]**2*R[2]**2)/(4*k*pi)


def test_example_1():
    """
    This is a sample example meant to test the basic functionality of
    electrostatics.py

    Consider a frame R, and a particle charge P(mass m, charge q) situated
    at its origin initially.
    Now, at time t=0, a field R.x+R.y+R.z is 'switched on' and kept that way
    for time 't0'.
    The user needs the x-coordinate of the particle charge as a function of
    time, after the field is switched off at time=t0.
    """
    
    #Basic initialization
    m, q = symbols('m q')
    P = ParticleCharge('P', m, q)
    P.set_motion(R, pos_vector = 0)
    field = R.x + R.y + R.z
    time = dynamicsymbols._t
    t0 = Symbol('t0')
    #The acceleration is equal to the electrostatic force experience by the
    #particle charge, divided by its mass
    acceleration = P.electrostatic_force(field)/P.mass
    #Use get_motion_acc from the mechanics core to find velocity and position
    #parameters using acceleration function and boundary conditions
    translation = get_motion_acc(acceleration, \
                                 0, P.pos_vector_wrt(R), frame = R)
    #Calculate the motion parameters of the particle charge at time t0
    acceleration = translation[0]
    velocity_at_t0 = translation[1].subs({time:t0})
    position_at_t0 = translation[2].subs({time:t0})
    #Set motion of P accordingly
    P.set_motion(R, trans_acc = 0, trans_vel_b = velocity_at_t0,
                 pos_vector_b = position_at_t0)
    #assert that the dot product of P's pos_vector wrt R, with R.x
    #matches the actual result
    assert P.pos_vector_wrt(R).dot(R.x) == \
           q*t*t0/m + q*t0**2/(2*m)


def test_example_2():
    """
    This is a sample example meant to test the basic functionality of
    electrostatics.py

    Consider a set of 3 particle charges of equal mass and charge placed
    symmetrically on a circle of unit radius around the origin, with one
    point lying at position vector of R.y
    We calculate the energy required to assemble this system of charges from
    infinity(with no kinetic energy in any particle) and find out values of
    the electrostatic fields and potentials at the origin
    """

    #Basic initialization
    from sympy import sin, cos, pi
    m, q = symbols('m q')
    P = ParticleCharge('P', m, q)
    Q = ParticleCharge('Q', m, q)
    S = ParticleCharge('S', m, q)
    #Place particle charges as per required configuration
    P.set_motion(R, pos_vector = R.y)
    Q.set_motion(R, pos_vector = -cos(pi/6) * R.x + -sin(pi/6) * R.y)
    S.set_motion(R, pos_vector = cos(pi/6) * R.x + -sin(pi/6) * R.y)
    #Check outputs of required values
    assert charge_assembly_energy(P, Q, S) == 3*3**(-0.5)*k*q**2
    assert S.electrostatic_field(R, 0) + \
           Q.electrostatic_field(R, 0) + \
           P.electrostatic_field(R, 0) == 0
    assert P.electrostatic_potential(R, 0) == \
           Q.electrostatic_potential(R, 0) == \
           S.electrostatic_potential(R, 0)
