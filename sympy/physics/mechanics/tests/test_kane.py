from __future__ import annotations
from sympy import (solve, lambdify, cos, expand, Matrix, sin, symbols, tan,
                   sqrt, S, zeros, eye)
from sympy.simplify.simplify import simplify
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                     RigidBody, KanesMethod, inertia, Particle,
                                     dot, find_dynamicsymbols)
from sympy.testing.pytest import raises, warns_deprecated_sympy, slow
from sympy.matrices.exceptions import NonInvertibleMatrixError


def test_invalid_coordinates():
    # Simple pendulum, but use symbols instead of dynamicsymbols
    l, m, g = symbols('l m g')
    q, u = symbols('q u')  # Generalized coordinate
    kd = [q.diff(dynamicsymbols._t) - u]
    N, O = ReferenceFrame('N'), Point('O')
    O.set_vel(N, 0)
    P = Particle('P', Point('P'), m)
    P.point.set_pos(O, l * (sin(q) * N.x - cos(q) * N.y))
    F = (P.point, -m * g * N.y)
    raises(ValueError, lambda: KanesMethod(N, [q], [u], kd, bodies=[P],
                                           forcelist=[F]))


def test_one_dof():
    # This is for a 1 dof spring-mass-damper case.
    # It is described in more detail in the KanesMethod docstring.
    q, u = dynamicsymbols('q u')
    qd, ud = dynamicsymbols('q u', 1)
    m, c, k = symbols('m c k')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, u * N.x)

    kd = [qd - u]
    FL = [(P, (-k * q - c * u) * N.x)]
    pa = Particle('pa', P, m)
    BL = [pa]

    KM = KanesMethod(N, [q], [u], kd)
    KM.kanes_equations(BL, FL)

    assert KM.frame == N
    assert KM.bodies == BL
    assert KM.loads == FL
    with warns_deprecated_sympy():
        assert KM.bodylist == BL
    with warns_deprecated_sympy():
        assert KM.forcelist == FL

    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand(-(q * k + u * c) / m)

    assert simplify(KM.rhs() -
                    KM.mass_matrix_full.LUsolve(KM.forcing_full)) == zeros(2, 1)

    assert (KM.linearize(A_and_B=True, )[0] == Matrix([[0, 1], [-k/m, -c/m]]))


def test_two_dof():
    # This is for a 2 d.o.f., 2 particle spring-mass-damper.
    # The first coordinate is the displacement of the first particle, and the
    # second is the relative displacement between the first and second
    # particles. Speeds are defined as the time derivatives of the particles.
    q1, q2, u1, u2 = dynamicsymbols('q1 q2 u1 u2')
    q1d, q2d, u1d, u2d = dynamicsymbols('q1 q2 u1 u2', 1)
    m, c1, c2, k1, k2 = symbols('m c1 c2 k1 k2')
    N = ReferenceFrame('N')
    P1 = Point('P1')
    P2 = Point('P2')
    P1.set_vel(N, u1 * N.x)
    P2.set_vel(N, (u1 + u2) * N.x)
    # Note we multiply the kinematic equation by an arbitrary factor
    # to test the implicit vs explicit kinematics attribute
    kd = [q1d/2 - u1/2, 2*q2d - 2*u2]

    # Now we create the list of forces, then assign properties to each
    # particle, then create a list of all particles.
    FL = [(P1, (-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2) * N.x), (P2, (-k2 *
        q2 - c2 * u2) * N.x)]
    pa1 = Particle('pa1', P1, m)
    pa2 = Particle('pa2', P2, m)
    BL = [pa1, pa2]

    # Finally we create the KanesMethod object, specify the inertial frame,
    # pass relevant information, and form Fr & Fr*. Then we calculate the mass
    # matrix and forcing terms, and finally solve for the udots.
    KM = KanesMethod(N, q_ind=[q1, q2], u_ind=[u1, u2], kd_eqs=kd)
    KM.kanes_equations(BL, FL)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand((-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2)/m)
    assert expand(rhs[1]) == expand((k1 * q1 + c1 * u1 - 2 * k2 * q2 - 2 *
                                    c2 * u2) / m)

    # Check that the explicit form is the default and kinematic mass matrix is identity
    assert KM.explicit_kinematics
    assert KM.mass_matrix_kin == eye(2)

    # Check that for the implicit form the mass matrix is not identity
    KM.explicit_kinematics = False
    assert KM.mass_matrix_kin == Matrix([[S(1)/2, 0], [0, 2]])

    # Check that whether using implicit or explicit kinematics the RHS
    # equations are consistent with the matrix form
    for explicit_kinematics in [False, True]:
        KM.explicit_kinematics = explicit_kinematics
        assert simplify(KM.rhs() -
                        KM.mass_matrix_full.LUsolve(KM.forcing_full)) == zeros(4, 1)

    # Make sure an error is raised if nonlinear kinematic differential
    # equations are supplied.
    kd = [q1d - u1**2, sin(q2d) - cos(u2)]
    raises(ValueError, lambda: KanesMethod(N, q_ind=[q1, q2],
                                           u_ind=[u1, u2], kd_eqs=kd))

def test_pend():
    q, u = dynamicsymbols('q u')
    qd, ud = dynamicsymbols('q u', 1)
    m, l, g = symbols('m l g')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, -l * u * sin(q) * N.x + l * u * cos(q) * N.y)
    kd = [qd - u]

    FL = [(P, m * g * N.x)]
    pa = Particle('pa', P, m)
    BL = [pa]

    KM = KanesMethod(N, [q], [u], kd)
    KM.kanes_equations(BL, FL)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    rhs.simplify()
    assert expand(rhs[0]) == expand(-g / l * sin(q))
    assert simplify(KM.rhs() -
                    KM.mass_matrix_full.LUsolve(KM.forcing_full)) == zeros(2, 1)


def test_rolling_disc():
    # Rolling Disc Example
    # Here the rolling disc is formed from the contact point up, removing the
    # need to introduce generalized speeds. Only 3 configuration and three
    # speed variables are need to describe this system, along with the disc's
    # mass and radius, and the local gravity (note that mass will drop out).
    q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1 q2 q3 u1 u2 u3')
    q1d, q2d, q3d, u1d, u2d, u3d = dynamicsymbols('q1 q2 q3 u1 u2 u3', 1)
    r, m, g = symbols('r m g')

    # The kinematics are formed by a series of simple rotations. Each simple
    # rotation creates a new frame, and the next rotation is defined by the new
    # frame's basis vectors. This example uses a 3-1-2 series of rotations, or
    # Z, X, Y series of rotations. Angular velocity for this is defined using
    # the second frame's basis (the lean frame).
    N = ReferenceFrame('N')
    Y = N.orientnew('Y', 'Axis', [q1, N.z])
    L = Y.orientnew('L', 'Axis', [q2, Y.x])
    R = L.orientnew('R', 'Axis', [q3, L.y])
    w_R_N_qd = R.ang_vel_in(N)
    R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)

    # This is the translational kinematics. We create a point with no velocity
    # in N; this is the contact point between the disc and ground. Next we form
    # the position vector from the contact point to the disc's center of mass.
    # Finally we form the velocity and acceleration of the disc.
    C = Point('C')
    C.set_vel(N, 0)
    Dmc = C.locatenew('Dmc', r * L.z)
    Dmc.v2pt_theory(C, N, R)

    # This is a simple way to form the inertia dyadic.
    I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)

    # Kinematic differential equations; how the generalized coordinate time
    # derivatives relate to generalized speeds.
    kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L]

    # Creation of the force list; it is the gravitational force at the mass
    # center of the disc. Then we create the disc by assigning a Point to the
    # center of mass attribute, a ReferenceFrame to the frame attribute, and mass
    # and inertia. Then we form the body list.
    ForceList = [(Dmc, - m * g * Y.z)]
    BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
    BodyList = [BodyD]

    # Finally we form the equations of motion, using the same steps we did
    # before. Specify inertial frame, supply generalized speeds, supply
    # kinematic differential equation dictionary, compute Fr from the force
    # list and Fr* from the body list, compute the mass matrix and forcing
    # terms, then solve for the u dots (time derivatives of the generalized
    # speeds).
    KM = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd)
    KM.kanes_equations(BodyList, ForceList)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    kdd = KM.kindiffdict()
    rhs = rhs.subs(kdd)
    rhs.simplify()
    assert rhs.expand() == Matrix([(6*u2*u3*r - u3**2*r*tan(q2) +
        4*g*sin(q2))/(5*r), -2*u1*u3/3, u1*(-2*u2 + u3*tan(q2))]).expand()
    assert simplify(KM.rhs() -
                    KM.mass_matrix_full.LUsolve(KM.forcing_full)) == zeros(6, 1)

    # This code tests our output vs. benchmark values. When r=g=m=1, the
    # critical speed (where all eigenvalues of the linearized equations are 0)
    # is 1 / sqrt(3) for the upright case.
    A = KM.linearize(A_and_B=True)[0]
    A_upright = A.subs({r: 1, g: 1, m: 1}).subs({q1: 0, q2: 0, q3: 0, u1: 0, u3: 0})
    import sympy
    assert sympy.sympify(A_upright.subs({u2: 1 / sqrt(3)})).eigenvals() == {S.Zero: 6}


def test_aux():
    # Same as above, except we have 2 auxiliary speeds for the ground contact
    # point, which is known to be zero. In one case, we go through then
    # substitute the aux. speeds in at the end (they are zero, as well as their
    # derivative), in the other case, we use the built-in auxiliary speed part
    # of KanesMethod. The equations from each should be the same.
    q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1 q2 q3 u1 u2 u3')
    q1d, q2d, q3d, u1d, u2d, u3d = dynamicsymbols('q1 q2 q3 u1 u2 u3', 1)
    u4, u5, f1, f2 = dynamicsymbols('u4, u5, f1, f2')
    u4d, u5d = dynamicsymbols('u4, u5', 1)
    r, m, g = symbols('r m g')

    N = ReferenceFrame('N')
    Y = N.orientnew('Y', 'Axis', [q1, N.z])
    L = Y.orientnew('L', 'Axis', [q2, Y.x])
    R = L.orientnew('R', 'Axis', [q3, L.y])
    w_R_N_qd = R.ang_vel_in(N)
    R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)

    C = Point('C')
    C.set_vel(N, u4 * L.x + u5 * (Y.z ^ L.x))
    Dmc = C.locatenew('Dmc', r * L.z)
    Dmc.v2pt_theory(C, N, R)
    Dmc.a2pt_theory(C, N, R)

    I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)

    kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L]

    ForceList = [(Dmc, - m * g * Y.z), (C, f1 * L.x + f2 * (Y.z ^ L.x))]
    BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
    BodyList = [BodyD]

    KM = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3, u4, u5],
                     kd_eqs=kd)
    (fr, frstar) = KM.kanes_equations(BodyList, ForceList)
    fr = fr.subs({u4d: 0, u5d: 0}).subs({u4: 0, u5: 0})
    frstar = frstar.subs({u4d: 0, u5d: 0}).subs({u4: 0, u5: 0})

    KM2 = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd,
                      u_auxiliary=[u4, u5])
    (fr2, frstar2) = KM2.kanes_equations(BodyList, ForceList)
    fr2 = fr2.subs({u4d: 0, u5d: 0}).subs({u4: 0, u5: 0})
    frstar2 = frstar2.subs({u4d: 0, u5d: 0}).subs({u4: 0, u5: 0})

    frstar.simplify()
    frstar2.simplify()

    assert (fr - fr2).expand() == Matrix([0, 0, 0, 0, 0])
    assert (frstar - frstar2).expand() == Matrix([0, 0, 0, 0, 0])


def test_parallel_axis():
    # This is for a 2 dof inverted pendulum on a cart.
    # This tests the parallel axis code in KanesMethod. The inertia of the
    # pendulum is defined about the hinge, not about the center of mass.

    # Defining the constants and knowns of the system
    gravity = symbols('g')
    k, ls = symbols('k ls')
    a, mA, mC = symbols('a mA mC')
    F = dynamicsymbols('F')
    Ix, Iy, Iz = symbols('Ix Iy Iz')

    # Declaring the Generalized coordinates and speeds
    q1, q2 = dynamicsymbols('q1 q2')
    q1d, q2d = dynamicsymbols('q1 q2', 1)
    u1, u2 = dynamicsymbols('u1 u2')
    u1d, u2d = dynamicsymbols('u1 u2', 1)

    # Creating reference frames
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')

    A.orient(N, 'Axis', [-q2, N.z])
    A.set_ang_vel(N, -u2 * N.z)

    # Origin of Newtonian reference frame
    O = Point('O')

    # Creating and Locating the positions of the cart, C, and the
    # center of mass of the pendulum, A
    C = O.locatenew('C', q1 * N.x)
    Ao = C.locatenew('Ao', a * A.y)

    # Defining velocities of the points
    O.set_vel(N, 0)
    C.set_vel(N, u1 * N.x)
    Ao.v2pt_theory(C, N, A)
    Cart = Particle('Cart', C, mC)
    Pendulum = RigidBody('Pendulum', Ao, A, mA, (inertia(A, Ix, Iy, Iz), C))

    # kinematical differential equations

    kindiffs = [q1d - u1, q2d - u2]

    bodyList = [Cart, Pendulum]

    forceList = [(Ao, -N.y * gravity * mA),
                 (C, -N.y * gravity * mC),
                 (C, -N.x * k * (q1 - ls)),
                 (C, N.x * F)]

    km = KanesMethod(N, [q1, q2], [u1, u2], kindiffs)
    (fr, frstar) = km.kanes_equations(bodyList, forceList)
    mm = km.mass_matrix_full
    assert mm[3, 3] == Iz

def test_input_format():
    # 1 dof problem from test_one_dof
    q, u = dynamicsymbols('q u')
    qd, ud = dynamicsymbols('q u', 1)
    m, c, k = symbols('m c k')
    N = ReferenceFrame('N')
    P = Point('P')
    P.set_vel(N, u * N.x)

    kd = [qd - u]
    FL = [(P, (-k * q - c * u) * N.x)]
    pa = Particle('pa', P, m)
    BL = [pa]

    KM = KanesMethod(N, [q], [u], kd)
    # test for input format kane.kanes_equations((body1, body2, particle1))
    assert KM.kanes_equations(BL)[0] == Matrix([0])
    # test for input format kane.kanes_equations(bodies=(body1, body 2), loads=(load1,load2))
    assert KM.kanes_equations(bodies=BL, loads=None)[0] == Matrix([0])
    # test for input format kane.kanes_equations(bodies=(body1, body 2), loads=None)
    assert KM.kanes_equations(BL, loads=None)[0] == Matrix([0])
    # test for input format kane.kanes_equations(bodies=(body1, body 2))
    assert KM.kanes_equations(BL)[0] == Matrix([0])
    # test for input format kane.kanes_equations(bodies=(body1, body2), loads=[])
    assert KM.kanes_equations(BL, [])[0] == Matrix([0])
    # test for error raised when a wrong force list (in this case a string) is provided
    raises(ValueError, lambda: KM._form_fr('bad input'))

    # 1 dof problem from test_one_dof with FL & BL in instance
    KM = KanesMethod(N, [q], [u], kd, bodies=BL, forcelist=FL)
    assert KM.kanes_equations()[0] == Matrix([-c*u - k*q])

    # 2 dof problem from test_two_dof
    q1, q2, u1, u2 = dynamicsymbols('q1 q2 u1 u2')
    q1d, q2d, u1d, u2d = dynamicsymbols('q1 q2 u1 u2', 1)
    m, c1, c2, k1, k2 = symbols('m c1 c2 k1 k2')
    N = ReferenceFrame('N')
    P1 = Point('P1')
    P2 = Point('P2')
    P1.set_vel(N, u1 * N.x)
    P2.set_vel(N, (u1 + u2) * N.x)
    kd = [q1d - u1, q2d - u2]

    FL = ((P1, (-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2) * N.x), (P2, (-k2 *
        q2 - c2 * u2) * N.x))
    pa1 = Particle('pa1', P1, m)
    pa2 = Particle('pa2', P2, m)
    BL = (pa1, pa2)

    KM = KanesMethod(N, q_ind=[q1, q2], u_ind=[u1, u2], kd_eqs=kd)
    # test for input format
    # kane.kanes_equations((body1, body2), (load1, load2))
    KM.kanes_equations(BL, FL)
    MM = KM.mass_matrix
    forcing = KM.forcing
    rhs = MM.inv() * forcing
    assert expand(rhs[0]) == expand((-k1 * q1 - c1 * u1 + k2 * q2 + c2 * u2)/m)
    assert expand(rhs[1]) == expand((k1 * q1 + c1 * u1 - 2 * k2 * q2 - 2 *
                                    c2 * u2) / m)


def test_implicit_kinematics():
    # Test that implicit kinematics can handle complicated equations that
    # explicit form struggles with See
    # https://github.com/sympy/sympy/issues/22626 This test was originally
    # designed to check whether the dynamical differential equations had
    # significantly fewer mathematical operations if the kinematical
    # differential equations were left in implicit form, but the reason the
    # equations had fewer operations was partially due to the fact that the
    # time differentiated holonomic constraint (quaternion length in this case)
    # was not provided to KanesMethod. At the time this test was written,
    # passing in a configuration constraint without passing in a velocity
    # constraint simply meant that the constraint was ignored and the
    # unconstrained dynamical differential equations were formed, which would
    # have much fewer operations than if the constraint was properly supplied.

    # Inertial frame
    NED = ReferenceFrame('NED')
    # Origin fixed in inertial reference frame
    NED_o = Point('NED_o')
    NED_o.set_vel(NED, 0)

    # Orientation generalized coordinates
    q_att = dynamicsymbols('lambda_0:4', real=True)  # quaternion scalars
    q_att_vec = Matrix(q_att)
    B = NED.orientnew('B', 'Quaternion', q_att)

    # Position generalized coordinates
    q_pos = dynamicsymbols('B_x:z')
    B_cm = NED_o.locatenew('B_cm', q_pos[0]*B.x + q_pos[1]*B.y + q_pos[2]*B.z)

    # Coordinates
    q_ind = q_att[1:] + q_pos  # lam1, lam2, lam3, B_x, B_y, B_z
    q_dep = [q_att[0]]  # lam0

    # Generalized speeds
    P, Q, R = dynamicsymbols('P, Q, R')
    U, V, W, X = dynamicsymbols('U, V,  W, X')

    u_ind = [U, V, W, P, Q, R]
    u_dep = [X]

    # Angular velocity
    B_ang_vel = B.ang_vel_in(NED)
    B.set_ang_vel(NED, P*B.x + Q*B.y + R*B.z)

    # Translational velocity
    B_cm_vel = B_cm.vel(NED)
    B_cm.set_vel(NED, U*B.x + V*B.y + W*B.z)

    # Kinematical differential equations
    # A generalized speed is needed for the fourth generalized coordinate to
    # comply with Kane's form of the definition of generalized speeds:
    # u = Y*q' where len(u) = len(q), and Y must be regular.
    # The quaternion kinematical differential equation takes this form:
    # [P] = Y*[lam0']
    # [Q]     [lam1']
    # [R]     [lam2']
    #         [lam3']
    # which is only 3 generalized speeds for 4 coordinates.
    # Solved it looks like:
    # dlam/dt = 1/2 * [0 -P -Q -R] * [lam0]
    #                 [P  0  R -Q]   [lam1]
    #                 [Q -R  0  P]   [lam2]
    #                 [R  Q -P  0]   [lam3]
    # which defines this kinematical differential equation of the four lam' as
    # a function of only 3 generalized speeds P, Q, R.
    B_ang_vel_kd = (B.ang_vel_in(NED) - B_ang_vel).simplify()
    B_ref_vel_kd = (B_cm.vel(NED) - B_cm_vel)

    # Add a differential equation that is equal to zero to manage the
    # additional generalized speed.
    kinematic_eqs = [
        (q_att_vec.T*q_att_vec.diff())[0, 0] - X,  # use velocity constraint
        B_ang_vel_kd.dot(B.x),
        B_ang_vel_kd.dot(B.y),
        B_ang_vel_kd.dot(B.z),
        B_ref_vel_kd.dot(B.x),
        B_ref_vel_kd.dot(B.y),
        B_ref_vel_kd.dot(B.z),
    ]
    qdots_sol = solve(kinematic_eqs, Matrix(q_dep + q_ind).diff())

    # Holonomic constraint: quaternion magnitude must equal 1
    # NOTE : You can analytically solve for lam0:
    # lam0 = sqrt(1 - lam1^2 - lam2^2 - lam3^2)
    # so the lam0 generalized coordinate could be eliminated.
    # [lam0, lam1, lam2, lam3]*[lam0] = lam0^2 + lam1^2 + lam2^2 + lam3^2 = 1
    #                          [lam1]
    #                          [lam2]
    #                          [lam3]
    config_cons = [(q_att_vec.T*q_att_vec)[0] - 1]  # unit norm
    # constraint is: lam0**2 + lam1**2 + lam2**2 + lam3**2 = 1
    # time derivative of the constraint is:
    # 2*lam0*lam0' + 2*lam1*lam1' + 2*lam2*lam2' + 2*lam3*lam3' = 0
    # or:
    # lam0*lam0' + lam1*lam1' + lam2*lam2' + lam3*lam3' = 0
    # [lam0, lam1, lam2, lam3]*[lam0'] = lam0*lam0' + ... + lam3*lam3' = X
    #                          [lam1']
    #                          [lam2']
    #                          [lam3']
    #vel_cons = [simplify((q_att_vec.T*q_att_vec.diff())[0, 0].xreplace(qdots_sol))]
    vel_cons = [X]

    # Single rigid body described by a central inertia dyadic
    M_B = symbols('M_B')
    J_B = inertia(B, *[S(f'J_B_{ax}')*(1 if ax[0] == ax[1] else -1)
                       for ax in ['xx', 'yy', 'zz', 'xy', 'yz', 'xz']])
    J_B = J_B.subs({S('J_B_xy'): 0, S('J_B_yz'): 0})
    RB = RigidBody('RB', B_cm, B, M_B, (J_B, B_cm))

    rigid_bodies = [RB]
    # Forces
    force_list = [
        # gravity pointing down
        (RB.masscenter, RB.mass*S('g')*NED.z),
        # generic forces and torques in body frame(inputs)
        (RB.frame, dynamicsymbols('T_z')*B.z),
        (RB.masscenter, dynamicsymbols('F_z')*B.z)
    ]

    # OPTION 1: Ignore the constraints and form the equations of motion as an
    # unconstrained system. The holonomic constraint would need to be handled
    # by augmenting with Lagrange multipliers.
    KM = KanesMethod(
        NED,
        q_dep + q_ind,
        u_dep + u_ind,
        kd_eqs=kinematic_eqs,
        bodies=rigid_bodies,
        forcelist=force_list,
    )
    KM.kanes_equations()

    # OPTION 2: Include the configuration constraint and let the velocity and
    # acceleration constraints be automatically handled.
    KM2 = KanesMethod(
        NED,
        q_ind,
        u_ind,
        q_dependent=q_dep,
        u_dependent=u_dep,
        configuration_constraints=config_cons,
        kd_eqs=kinematic_eqs,
        bodies=rigid_bodies,
        forcelist=force_list,
    )
    KM2.kanes_equations()

    # OPTION 3: Provide both the configuration and (simplified) velocity
    # constraints.
    KM3 = KanesMethod(
        NED,
        q_ind,
        u_ind,
        q_dependent=q_dep,
        u_dependent=u_dep,
        configuration_constraints=config_cons,
        velocity_constraints=vel_cons,
        kd_eqs=kinematic_eqs,
        bodies=rigid_bodies,
        forcelist=force_list,
    )
    KM3.kanes_equations()

    # TODO : OPTION 4: Solve for the dependent coordinate analytically.

    assert simplify(KM2.forcing - KM3.forcing) == zeros(len(KM2.u), 1)
    assert KM2.mass_matrix[:-1, :] == KM3.mass_matrix[:-1, :]
    assert KM2.mass_matrix[-1, :] == 2*KM3.mass_matrix[-1, :]
    assert KM2.forcing[-1] == 0

    n_ops_explicit = sum(
        [x.count_ops() for x in KM.forcing_full] +
        [x.count_ops() for x in KM.mass_matrix_full]
    )
    #mass_matrix_kin_explicit = KM.mass_matrix_kin
    forcing_kin_explicit = KM.forcing_kin

    # The effect of changing this attribute to False is only that Mk*q' = Fk is
    # combined with the dynamical differential equations instead of I*q' =
    # Inv(Mk)*Fk, making the kinematical differential equations shorter.
    KM.explicit_kinematics = False
    n_ops_implicit = sum(
        [x.count_ops() for x in KM.forcing_full] +
        [x.count_ops() for x in KM.mass_matrix_full]
    )
    # Save implicit kinematic matrices to use later
    mass_matrix_kin_implicit = KM.mass_matrix_kin
    forcing_kin_implicit = KM.forcing_kin

    # Expecting implicit form to be less than 15% of the flops
    assert n_ops_implicit / n_ops_explicit < .02

    # Ideally we would check that implicit and explicit equations give the same
    # result as done in test_one_dof, but the whole raison-d'etre of the
    # implicit equations is to deal with problems such as this one where the
    # explicit form is too complicated to handle, especially the angular part
    # (i.e. tests would be too slow)
    # Instead, we check that the kinematic equations are correct using more
    # fundamental tests:
    #
    # (1) that we recover the kinematic equations we have provided
    # M_imp*q' = F_imp
    expected = Matrix(kinematic_eqs)
    calculated = mass_matrix_kin_implicit*KM.q.diff() - forcing_kin_implicit
    assert simplify(expected - calculated) == zeros(len(KM.q), 1)

    # M_exp*q' = F_exp
    expected = Matrix([qdots_sol[qdi] for qdi in KM.q.diff()])
    calculated = forcing_kin_explicit
    # TODO : Slow comparison.
    #assert simplify(expected - calculated) == zeros(len(KM.q), 1)
    # Alternative numerical check:
    syms = find_dynamicsymbols(expected - calculated)
    eval_comp = lambdify(list(syms), expected - calculated)
    for v in eval_comp(*list(range(1, len(syms) + 1)))[:]:
        if abs(v) > 1e-13:
            raise AssertionError

    # (2) that rate of the quaternion matches the textbook solution
    # NOTE : The diagonal is not zero because we introduce a generalized speed
    # X = vel_con = 0.
    qdot_candidate = Matrix([qdots_sol[qdi] for qdi in q_att_vec.diff()])
    qdot_candidate = qdot_candidate.subs(2*config_cons[0] + 1, 1)

    quat_dot_textbook = Matrix([
        [2*X, -P, -Q, -R],
        [P,  2*X,  R, -Q],
        [Q, -R,  2*X,  P],
        [R,  Q, -P,  2*X],
    ])*q_att_vec/2

    assert simplify(qdot_candidate - quat_dot_textbook) == zeros(4, 1)

    # sub the config constraint in the candidate solution and compare to the
    # implicit rhs
    lambda_0_sol = solve(config_cons[0], q_att_vec[0])[1]
    # This seems to do M_imp*F_exp == F_imp, why???
    lhs_candidate = simplify(
        (mass_matrix_kin_implicit[0:4, 0:4]*qdot_candidate).subs(
            {q_att_vec[0]: lambda_0_sol}))
    assert simplify(
        lhs_candidate - forcing_kin_implicit[[0, 1, 2, 3], :]
    ) == zeros(4, 1)


@slow
def test_issue_24887():
    # Spherical pendulum
    g, l, m, c = symbols('g l m c')
    q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1:4 u1:4')
    N = ReferenceFrame('N')
    A = ReferenceFrame('A')
    A.orient_body_fixed(N, (q1, q2, q3), 'zxy')
    N_w_A = A.ang_vel_in(N)
    # A.set_ang_vel(N, u1 * A.x + u2 * A.y + u3 * A.z)
    kdes = [N_w_A.dot(A.x) - u1, N_w_A.dot(A.y) - u2, N_w_A.dot(A.z) - u3]
    O = Point('O')
    O.set_vel(N, 0)
    Po = O.locatenew('Po', -l * A.y)
    Po.set_vel(A, 0)
    P = Particle('P', Po, m)
    kane = KanesMethod(N, [q1, q2, q3], [u1, u2, u3], kdes, bodies=[P],
                       forcelist=[(Po, -m * g * N.y)])
    kane.kanes_equations()
    expected_md = m * l ** 2 * Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 1]])
    expected_fd = Matrix([
        [l*m*(g*(sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3)) - l*u2*u3)],
        [0], [l*m*(-g*(sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1)) + l*u1*u2)]])
    assert find_dynamicsymbols(kane.forcing).issubset({q1, q2, q3, u1, u2, u3})
    assert simplify(kane.mass_matrix - expected_md) == zeros(3, 3)
    assert simplify(kane.forcing - expected_fd) == zeros(3, 1)

def test_Ars_warning():
    """
    Tests that the additional warning is raised when the mass matrix is
    singular and that the correct error message is displayed.
    It is a simple pendulum on a post, where the pole must not move in the
    y direction. Clearly, if y is not the dependent speed,
    Ars will be singular.
    """
    N, A = symbols('N A', cls=ReferenceFrame)
    O, P0, P1, P2 = symbols('O P0 P1 P2', cls=Point)
    t = dynamicsymbols._t
    O.set_vel(N, 0)

    l = symbols('l')

    x, y, q = dynamicsymbols('x y q')
    ux, uy, uq = dynamicsymbols('ux uy uq')

    A.orient_axis(N, q, N.z)

    P0.set_pos(O, x * N.x + y * N.y)
    P0.set_vel(N, ux * N.x + uy * N.y)
    P1.set_pos(P0, 2*l * N.y)
    P2.set_pos(P1, l * A.x)
    P2.v2pt_theory(P1, N, A)

    kd = Matrix([ux - x.diff(), uy - y.diff(), uq - q.diff()])

    config_constr = Matrix([y])
    speed_constr = config_constr.diff(t)

    q_ind = [x, y]
    q_dep = [q]
    u_ind = [ux, uy]
    u_dep = [uq]

    msg = ("Matrix det == 0; not invertible. Ars is singular, maybe due to a "
           "wrong set of dependent speeds. Try a different set of "
           "dependent speeds.")

    with raises(NonInvertibleMatrixError) as excinfo:
        _ = KanesMethod(N,
                           q_ind,
                           u_ind,
                           kd_eqs=kd,
                           q_dependent=q_dep,
                           u_dependent=u_dep,
                           configuration_constraints=config_constr,
                           velocity_constraints=speed_constr
        )

    assert msg in str(excinfo.value)

def test_constraint_combos():
    # This tests the combinations of ways you can supply the constraints.

    l, m = symbols('l m')
    q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1:4 u1:4')
    N, A = ReferenceFrame('N'), ReferenceFrame('A')
    A.orient_body_fixed(N, (q1, q2, q3), 'zxy')
    N_w_A = A.ang_vel_in(N)
    kdes = [N_w_A.dot(A.x) - u1, N_w_A.dot(A.y) - u2, N_w_A.dot(A.z) - u3]
    O = Point('O')
    O.set_vel(N, 0)
    Po = O.locatenew('Po', -l*A.y)
    Po.set_vel(A, 0)
    P = Particle('P', Po, m)

    # no constraints
    kane = KanesMethod(
        N,
        [q1, q2, q3],
        [u1, u2, u3],
        kdes, bodies=[P]
    )
    assert kane.holonomic_constraints == Matrix([])
    assert kane.nonholonomic_constraints == Matrix([])
    assert kane.velocity_constraints == Matrix([])

    # empty constraints
    kane = KanesMethod(
        N,
        [q1, q2, q3],
        [u1, u2, u3],
        kdes, bodies=[P],
        configuration_constraints=[],
        nonholonomic_constraints=[],
    )
    assert kane.holonomic_constraints == Matrix([])
    assert kane.nonholonomic_constraints == Matrix([])
    assert kane.velocity_constraints == Matrix([])

    # only holonomic constraints
    kane = KanesMethod(
        N,
        [q2, q3],
        [u2, u3],
        kdes, bodies=[P],
        q_dependent=[q1],
        configuration_constraints=[q1 - q2],
        u_dependent=[u1],
    )
    assert kane.holonomic_constraints == Matrix([q1 - q2])
    assert kane.nonholonomic_constraints == Matrix([])
    assert kane.velocity_constraints == Matrix([kane.kindiffdict()[q1.diff()] -
                                                kane.kindiffdict()[q2.diff()]])

    # only nonholonomic constraints
    kane = KanesMethod(
        N,
        [q1, q2, q3],
        [u2, u3],
        kdes, bodies=[P],
        nonholonomic_constraints=[u1 - u2],
        u_dependent=[u1],
    )
    assert kane.holonomic_constraints == Matrix([])
    assert kane.nonholonomic_constraints == Matrix([u1 - u2])
    assert kane.velocity_constraints == Matrix([u1 - u2])

    # holonomic and nonholonomic constraints
    kane = KanesMethod(
        N,
        [q2, q3],
        [u2],
        kdes, bodies=[P],
        configuration_constraints=[q1 - q2],
        q_dependent=[q1],
        nonholonomic_constraints=[u2 - u3],
        u_dependent=[u1, u3],
    )
    assert kane.holonomic_constraints == Matrix([q1 - q2])
    assert kane.nonholonomic_constraints == Matrix([u2 - u3])
    assert kane.velocity_constraints == Matrix([kane.kindiffdict()[q1.diff()] -
                                                kane.kindiffdict()[q2.diff()],
                                                u2 - u3])

    # holonomic and velocity constraints
    kane = KanesMethod(
        N,
        [q2, q3],
        [u2],
        kdes, bodies=[P],
        configuration_constraints=[q1 - q2],
        q_dependent=[q1],
        velocity_constraints=[q1.diff() - q2.diff(), u2 - u3],
        u_dependent=[u1, u3],
    )
    assert kane.holonomic_constraints == Matrix([q1 - q2])
    assert kane.nonholonomic_constraints == Matrix([])
    assert kane.velocity_constraints == Matrix([kane.kindiffdict()[q1.diff()] -
                                                kane.kindiffdict()[q2.diff()],
                                                u2 - u3])

    # nonholonomic and velocity not simultaneously allowed
    with raises(ValueError):
        kane = KanesMethod(
            N,
            [q1, q2, q3],
            [u1, u2, u3],
            kdes, bodies=[P],
            velocity_constraints=[],
            nonholonomic_constraints=[],
        )

    # incorrect number of dependent coordinates
    with raises(ValueError):
        kane = KanesMethod(
            N,
            [q1, q2, q3],
            [u2],
            kdes, bodies=[P],
            configuration_constraints=[q1 - q2],
            q_dependent=[],
            nonholonomic_constraints=[u2 - u3],
            u_dependent=[u1, u3],
        )

    # incorrect number of dependent speeds
    with raises(ValueError):
        kane = KanesMethod(
            N,
            [q1, q2, q3],
            [u1, u2],
            kdes, bodies=[P],
            configuration_constraints=[q1 - q2],
            q_dependent=[q1],
            nonholonomic_constraints=[u2 - u3],
            u_dependent=[u3],
        )
