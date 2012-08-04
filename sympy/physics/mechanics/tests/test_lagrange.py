from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                    RigidBody, LagrangesMethod, Kane, Particle)
from sympy import symbols, pi, sin, cos, simplify

def test_disc_on_an_incline_plane():
    # Disc rolling on an inclined plane
    # First the generalized coordinates are created. The mass center of the
    # disc is located from top vertex of the inclined plane by the generalized
    # coordinate 'y'. The orientation of the disc is defined by the angle
    # 'theta'. The mass of the disc is 'm' and its radius is 'R'. The length of
    # the inclined path is 'l', the angle of inclination is 'alpha'. 'g' is the
    # gravitational constant.
    y, theta = dynamicsymbols('y theta')
    yd, thetad = dynamicsymbols('y theta', 1)
    m, g, R, l, alpha = symbols('m g R l alpha')

    # Next, we create the inertial reference frame 'N'. A reference frame 'A'
    # is attached to the inclined plane. Finally a frame is created which is attached to the disk.
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [pi/2 -alpha, N.z])
    B = A.orientnew('B', 'Axis', [-theta, A.z])

    # Creating the disc 'D'; we create the point that represents the mass
    # center of the disc and set its velocity. The inertia dyadic of the disc
    # is created. Finally, we create the disc.
    Do = Point('Do')
    Do.set_vel(N, yd * A.x)
    I = m * R**2 / 2 * B.z | B.z
    D = RigidBody('D', Do, B, m, (I, Do))

    # To construct the Lagrangian, 'L', of the disc, we determine its kinetic
    # and potential energies, T and U, respectively. L is defined as the
    # difference between T and U.
    T = (m *(Do.vel(N) & Do.vel(N)))/2 + (B.ang_vel_in(N) & (I & B.ang_vel_in(N)))/2
    U = m * g * (l -y) * sin(alpha)
    L = T - U

    # We then create the list of generalized coordinates and constraint
    # equations. The constraint arises due to the disc rolling without slip on
    # on the inclined path. Also, the constraint is holonomic but we supply the
    # differentiated holonomic equation as the 'LagrangesMethod' class requires
    # that. We then invoke the 'LagrangesMethod' class and supply it the
    # necessary arguments and generate the equations of motion. The'rhs' method
    # solves for the q_double_dots (i.e. the second derivative with respect to
    # time  of the generalized coordinates and the lagrange multiplers.
    q = [y, theta]
    coneq = [yd - R * thetad]
    m = LagrangesMethod(L, q, coneq)
    m.form_lagranges_equations()
    assert m.rhs("GE")[2] ==  2*g*sin(alpha)/3

def test_simp_pen():
    # This test compares that the equations from Kane and LagrangesMethod are
    # indeed equivalent. The system under consideration is the simple pendulum.
    # The first segment of the code is for the 'Kane' class. The equations are
    # derived by the 'LagrangesMethod' class following that.
    # We begin by creating the generalized coordinates and generalized speeds
    # as per the requirements of 'Kane'. Also we created the associate symbols
    # that characterize the system: 'm' is the mass of the bob, l is the length
    # of the massless string connecting the bob to a point O fixed in the
    # inertial frame.
    q, u = dynamicsymbols('q u')
    qd, ud = dynamicsymbols('q u ', 1)
    l, m, g = symbols('l m g')

    # We then create the inertial frame and a frame attached to the massless
    # string following which we define the inertial angular velocity of the
    # string.
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [q, N.z])
    A.set_ang_vel(N, u * N.z)

    # Next, we create the point O and fix it in the inertial frame. We then
    # locate the point P to which the bob is attached. Its corresponding
    # velocity is then determined by the 'two point formula'.
    O = Point('O')
    O.set_vel(N, 0)
    P = O.locatenew('P', l * A.x)
    P.v2pt_theory(O, N, A)

    # The 'Particle' which represents the bob is then created.
    Pa = Particle('Pa', P, m)

    # We go through the formalism as described by the 'Kane' class to generate
    # the equations of motion. Please refer to the documentation for
    # clarification on the procedure here.
    kd = [qd - u]
    FL = [(P, m * g * N.x)]
    BL = [Pa]

    KM = Kane(N)
    KM.coords([q])
    KM.speeds([u])
    KM.kindiffeq(kd)

    (fr, frstar) = KM.kanes_equations(FL, BL)
    mm = KM.mass_matrix_full
    fo = KM.forcing_full
    qudots = mm.inv() * fo

    # This section is pertinent to the 'LagrangesMethod' class. The orientation
    # of the string fixed frame is redefined as 'LagrangesMethod' doesn't require
    # generalized speeds, per se. (Lagrangian mechanics requires 'simple'
    # generalized speeds)
    A.orientnew('A', 'Axis', [q, N.z])
    A.set_ang_vel(N, qd *A.z)
    P.v2pt_theory(O,N,A)
    T = 1/2.0 * m * P.vel(N) & P.vel(N) # T is the kinetic energy of the system
    V = - m * g * l * cos(q) # V is the potential energy of the system
    L = T - V # L is the Lagrangian

    # The 'LagrangesMethod' class is invoked and the equations of motion are generated.
    l = LagrangesMethod(L, [q])
    l.form_lagranges_equations()
    l.rhs("GE")

    # Finally, we the results of both methods is compared and it's seen that
    # they are indeed equivalent.
    assert l.rhs("GE")[0] == qudots.subs({u : qd})[0]
    assert l.rhs("GE")[1] == qudots.subs({u : qd})[1]

def test_dub_pen():

    # Yet another test to compare that the equations generated by LagrangesMethod
    # are equivalent. The system considered is the double pendulum. Like in the
    # test of the simple pendulum above, we begin by creating the generalized
    # coordinates and generalized speeds (associated with 'Kane'). Also created are
    # the simple generalized speeds and accelerations which will be used later.
    # Following this we create frames and points necessary for the kinematics.
    # The procedure isn't explicitly explained as this is similar to the simple
    # pendulum. Also this is documented on the pydy.org website.

    q1, q2 = dynamicsymbols('q1 q2')
    q1d, q2d = dynamicsymbols('q1 q2', 1)
    q1dd, q2dd = dynamicsymbols('q1 q2', 2)
    u1, u2 = dynamicsymbols('u1 u2')
    u1d, u2d = dynamicsymbols('u1 u2', 1)
    l, m, g = symbols('l m g')

    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [q1, N.z])
    B = N.orientnew('B', 'Axis', [q2, N.z])

    A.set_ang_vel(N, u1 * N.z)
    B.set_ang_vel(N, u2 * N.z)

    O = Point('O')
    P = O.locatenew('P', l * A.x)
    R = P.locatenew('R', l * B.x)

    O.set_vel(N, 0)
    P.v2pt_theory(O, N, A)
    R.v2pt_theory(P, N, B)

    ParP = Particle('ParP', P, m)
    ParR = Particle('ParR', R, m)

    kd = [q1d - u1, q2d - u2]
    FL = [(P, m * g * N.x), (R, m * g * N.x)]
    BL = [ParP, ParR]

    # Here we go through the steps to generate Kane's equations.
    KM = Kane(N)
    KM.coords([q1, q2])
    KM.speeds([u1, u2])
    KM.kindiffeq(kd)
    (fr, frstar) = KM.kanes_equations(FL, BL)

    # Here we go through the process to generate equations by LagrangesMethod
    A.set_ang_vel(N, q1d * A.z)
    B.set_ang_vel(N, q2d * A.z)
    P.v2pt_theory(O, N, A)
    R.v2pt_theory(P, N, B)

    # The kinetic energy of the system
    T = (m * P.vel(N) & P.vel(N)/2.0) + (m * R.vel(N) & R.vel(N)/2.0)
    # V1 and V2 are the potential energies of particles 'ParP' and 'ParR'
    V1 = -m * g * l * cos(q1)
    V2 = (- m * g * l * cos(q1) - m * g * l * cos(q2))
    V = V1 + V2
    # Then the Lagrangian is given by
    L = T - V
    l = LagrangesMethod(L, [q1, q2])
    l.form_lagranges_equations()

    # Checking that the 'RHS' generated by both methods are identical.
    # Because there were some issues with how 'Derivative' works with
    # 'dynamicsymbols', the 'dynamicsymbols' objects are replaced by 'symbols'
    # in the RHS vectors.

    Q1d, Q2d, Q1dd, Q2dd = symbols('Q1d Q2d Q1dd Q2dd')
    qd_dict = {q1d : Q1d, q2d : Q2d}
    qdd_dict = {q1dd : Q1dd, q2dd : Q2dd}
    u_dict = {u1 : Q1d, u2 : Q2d}
    ud_dict = {u1d : Q1dd, u2d : Q2dd}

    KMrhs = KM.mass_matrix.inv() * KM.forcing
    LMrhs = l.mass_matrix.inv() * l.forcing

    comp_mat = KMrhs - LMrhs
    comp_mat_sub = ((((comp_mat.subs(qd_dict)).subs(qdd_dict).subs(u_dict).subs(ud_dict))))
    assert simplify(comp_mat_sub[0]) == 0
    assert simplify(comp_mat_sub[0]) == 0
