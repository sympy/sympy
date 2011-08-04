from sympy import (symbols, Matrix, sin, cos, tan, expand,
                   solve_linear_system_LU)
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                     RigidBody, Kane, inertia, Particle)

def test_one_dof():
    # This is for a 1 dof spring-mass-damper case.
    # It is described in more detail in the Kane docstring.
    q, qd, u, ud = dynamicsymbols('q qd u ud')
    m, c, k = symbols('m c k')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, u * N.x)

    kd = [qd - u]
    FL = [(P, (-k * q - c * u) * N.x)]
    pa = Particle()
    pa.mass = m
    pa.point = P
    BL = [pa]

    KM = Kane(N)
    KM.coords([q])
    KM.speeds([u])
    KM.kindiffeq(kd)
    KM.form_fr(FL)
    KM.form_frstar(BL)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand(-(q * k + u * c) / m)

def test_two_dof():
    # This is for a 2 d.o.f., 2 particle spring-mass-damper.
    # The first coordinate is the displacement of the first particle, and the
    # second is the relative displacement between the first and second
    # particles. Speeds are defined as the time derivatives of the particles.
    q1, q2, q1d, q2d = dynamicsymbols('q1 q2 q1d q2d')
    u1, u2, u1d, u2d = dynamicsymbols('u1 u2 u1d u2d')
    m, c1, c2, k1, k2 = symbols('m c1 c2 k1 k2')
    N = ReferenceFrame('N')
    P1 = Point('P1')
    P2 = Point('P2')
    P1.set_vel(N, u1 * N.x)
    P2.set_vel(N, (u1 + u2) * N.x)
    kd = [q1d - u1, q2d - u2]

    # Now we create the list of forces, then assign properties to each
    # particle, then create a list of all particles.
    FL = [(P1, (-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2) * N.x), (P2, (-k2 *
        q2 - c2 * u2) * N.x)]
    pa1 = Particle()
    pa2 = Particle()
    pa1.mass = m
    pa2.mass = m
    pa1.point = P1
    pa2.point = P2
    BL = [pa1, pa2]

    # Finally we create the Kane object, specify the inertial frame, pass
    # relevant information, and form Fr & Fr*. Then we calculate the mass
    # matrix and forcing terms, and finally solve for the udots.
    KM = Kane(N)
    KM.coords([q1, q2])
    KM.speeds([u1, u2])
    KM.kindiffeq(kd)
    KM.form_fr(FL)
    KM.form_frstar(BL)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand((-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2)/m)
    assert expand(rhs[1]) == expand((k1 * q1 + c1 * u1 - 2 * k2 * q2 - 2 *
                                    c2 * u2) / m)

def test_pend():
    q, qd, u, ud = dynamicsymbols('q qd u ud')
    m, l, g = symbols('m l g')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, -l * u * sin(q) * N.x + l * u * cos(q) * N.y)
    kd = [qd - u]

    FL = [(P, m * g * N.x)]
    pa = Particle()
    pa.mass = m
    pa.point = P
    BL = [pa]

    KM = Kane(N)
    KM.coords([q])
    KM.speeds([u])
    KM.kindiffeq(kd)
    KM.form_fr(FL)
    KM.form_frstar(BL)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    rhs.simplify()
    assert expand(rhs[0]) == expand(-g / l * sin(q))


def test_rolling_disc():
    # Rolling Disc Example
    # Here the rolling disc is formed from the contact point up, removing the
    # need to introduce generalized speeds. Only 3 configuration and three
    # speed variables are need to describe this system, along with the disc's
    # mass and radius, and the local graivty (note that mass will drop out).
    q1, q2, q3, q1d, q2d, q3d = dynamicsymbols('q1 q2 q3 q1d q2d q3d')
    u1, u2, u3, u1d, u2d, u3d = dynamicsymbols('u1 u2 u3 u1d u2d u3d')
    r, m, g = symbols('r m g')

    # The kinematics are formed by a series of simple rotations. Each simple
    # rotation creates a new frame, and the next rotation is defined by the new
    # frame's basis vectors. This example uses a 3-1-2 series of rotations, or
    # Z, X, Y series of rotations. Angular velocity for this is defined using
    # the second frame's basis (the lean frame).
    N = ReferenceFrame('N')
    Y = N.orientnew('Y', 'Simple', q1, 3)
    L = Y.orientnew('L', 'Simple', q2, 1)
    R = L.orientnew('R', 'Simple', q3, 2)
    R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
    R.set_ang_acc(N, R.ang_vel_in(N).dt(R) + (R.ang_vel_in(N) ^ R.ang_vel_in(N)))

    # This is the translational kinematics. We create a point with no velocity
    # in N; this is the contact point between the disc and ground. Next we form
    # the position vector from the contact point to the disc mass center.
    # Finally we form the velocity and acceleration of the disc.
    C = Point('C')
    C.set_vel(N, 0)
    Dmc = C.locatenew('Dmc', r * L.z)
    Dmc.v2pt_theory(C, N, R)
    Dmc.a2pt_theory(C, N, R)

    # This is a simple way to form the inertia dyadic.
    I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)

    # Kinematic differential equations; how the generalized coordinate time
    # derivatives relate to generalized speeds.
    kd = [q1d - u3/cos(q3), q2d - u1, q3d - u2 + u3 * tan(q2)]

    # Creation of the force list; it is the gravitational force at the mass
    # center of the disc. Then we create the disc by assigning a Point to the
    # mass center attribute, a ReferenceFrame to the frame attribute, and mass
    # and inertia. Then we form the body list.
    ForceList = [(Dmc, - m * g * Y.z)]
    BodyD = RigidBody()
    BodyD.mc = Dmc
    BodyD.inertia = (I, Dmc)
    BodyD.frame = R
    BodyD.mass = m
    BodyList = [BodyD]

    # Finally we form the equations of motion, using the same steps we did
    # before. Specify inertial frame, supply generalized speeds, supply
    # kinematic differential equation dictionary, compute Fr from the force
    # list and Fr* fromt the body list, compute the mass matrix and forcing
    # terms, then solve for the u dots (time derivatives of the generalized
    # speeds).
    KM = Kane(N)
    KM.coords([q1, q2, q3])
    KM.speeds([u1, u2, u3])
    KM.kindiffeq(kd)
    fr = KM.form_fr(ForceList)
    frstar = KM.form_frstar(BodyList)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    sub_dict = solve_linear_system_LU(Matrix([KM._k_kqdot.T,
        -(KM._k_ku*Matrix(KM._u) + KM._f_k).T]).T, KM._qdot)
    rhs = rhs.subs(sub_dict)
    assert rhs.expand() == Matrix([(10*u2*u3*r - 5*u3**2*r*tan(q2) +
        4*g*sin(q2))/(5*r), -2*u1*u3/3, u1*(-2*u2 + u3*tan(q2))]).expand()

