from sympy import symbols, Matrix, solve, simplify, cos, sin, atan, sqrt
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point,\
        dot, cross, inertia, KanesMethod, Particle, RigidBody

from sympy.physics.mechanics.linearize import Linearizer

def test_linearize_rolling_disc():
    # Symbols for time and constant parameters
    t, r, m, g, v = symbols('t r m g v')

    # Configuration variables and their time derivatives
    q1, q2, q3, q4, q5, q6 = q = dynamicsymbols('q1:7')
    q1d, q2d, q3d, q4d, q5d, q6d = qd = [qi.diff(t) for qi in q]

    # Generalized speeds and their time derivatives
    u = dynamicsymbols('u:6')
    u1, u2, u3, u4, u5, u6 = u = dynamicsymbols('u1:7')
    u1d, u2d, u3d, u4d, u5d, u6d = ud = [ui.diff(t) for ui in u]

    # Reference frames
    N = ReferenceFrame('N')                   # Inertial frame
    NO = Point('NO')                          # Inertial origin
    A = N.orientnew('A', 'Axis', [q1, N.z])   # Yaw intermediate frame
    B = A.orientnew('B', 'Axis', [q2, A.x])   # Lean intermediate frame
    C = B.orientnew('C', 'Axis', [q3, B.y])   # Disc fixed frame
    CO = NO.locatenew('CO', q4*N.x + q5*N.y + q6*N.z)      # Disc center

    # Disc angular velocity in N expressed using time derivatives of coordinates
    w_c_n_qd = C.ang_vel_in(N)
    w_b_n_qd = B.ang_vel_in(N)

    # Inertial angular velocity and angular acceleration of disc fixed frame
    C.set_ang_vel(N, u1*B.x + u2*B.y + u3*B.z)

    # Disc center velocity in N expressed using time derivatives of coordinates
    v_co_n_qd = CO.pos_from(NO).dt(N)

    # Disc center velocity in N expressed using generalized speeds
    CO.set_vel(N, u4*C.x + u5*C.y + u6*C.z)

    # Disc Ground Contact Point
    P = CO.locatenew('P', r*B.z)
    P.v2pt_theory(CO, N, C)

    # Configuration constraint
    f_c = Matrix([q6 - dot(CO.pos_from(P), N.z)])

    # Velocity level constraints
    f_v = Matrix([dot(P.vel(N), uv) for uv in C])

    # Kinematic differential equations
    kindiffs = Matrix([dot(w_c_n_qd - C.ang_vel_in(N), uv) for uv in B] +
                        [dot(v_co_n_qd - CO.vel(N), uv) for uv in N])
    qdots = solve(kindiffs, qd)

    # Set angular velocity of remaining frames
    B.set_ang_vel(N, w_b_n_qd.subs(qdots))
    C.set_ang_acc(N, C.ang_vel_in(N).dt(B) + cross(B.ang_vel_in(N), C.ang_vel_in(N)))

    # Active forces
    F_CO = m*g*A.z

    # Inertia force
    R_star_CO = -m*CO.acc(N)

    # Create inertia dyadic of disc C about point CO
    I = (m * r**2) / 4
    J = (m * r**2) / 2
    I_C_CO = inertia(C, I, J, I)

    Disc = RigidBody('Disc', CO, C, m, (I_C_CO, CO))
    BL = [Disc]
    FL = [(CO, F_CO)]
    KM = KanesMethod(N, [q1, q2, q3, q4, q5], [u1, u2, u3], kd_eqs=kindiffs,
            q_dependent=[q6], configuration_constraints=f_c,
            u_dependent=[u4, u5, u6], velocity_constraints=f_v)
    (fr, fr_star) = KM.kanes_equations(FL, BL)

    # Test generalized form equations
    gen_form = KM._general_form()
    assert gen_form[0] == f_c
    assert gen_form[1] == f_v
    assert gen_form[2] == f_v.diff(t)
    sol = solve(gen_form[3] + gen_form[4], qd)
    for qi in qd:
        assert sol[qi] == qdots[qi]
    assert simplify(gen_form[5] + gen_form[6] - fr - fr_star) == Matrix([0, 0, 0])

    # Perform the linearization
    # Precomputed trim condition:
    eq_q = {q6: -r*cos(q2)}
    eq_u = {u1: 0,
            u2: sin(q2)*q1d + q3d,
            u3: cos(q2)*q1d,
            u4: -r*(sin(q2)*q1d + q3d)*cos(q3),
            u5: 0,
            u6: -r*(sin(q2)*q1d + q3d)*sin(q3)}
    eq_qd = {q2d: 0,
            q4d: -r*(sin(q2)*q1d + q3d)*cos(q1),
            q5d: -r*(sin(q2)*q1d + q3d)*sin(q1),
            q6d: 0}
    eq_ud = {u1d: 4*g*sin(q2)/(5*r) + sin(2*q2)*q1d**2/2 + 6*cos(q2)*q1d*q3d/5,
            u2d: 0,
            u3d: 0,
            u4d: r*(sin(q2)*sin(q3)*q1d*q3d + sin(q3)*q3d**2),
            u5d: r*(4*g*sin(q2)/(5*r) + sin(2*q2)*q1d**2/2 + 6*cos(q2)*q1d*q3d/5),
            u6d: -r*(sin(q2)*cos(q3)*q1d*q3d + cos(q3)*q3d**2)}

    linearizer = KM.to_linearizer()
    A, B = linearizer.linearize(eq_q, eq_u, eq_qd, eq_ud, A_and_B=True)

    upright_nominal = {q1d: 0, q2: 0, m: 1, r: 1, g: 1}

    # Precomputed solution
    A_sol = Matrix([[0, 0, 0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 1, 0],
                    [sin(q1)*q3d, 0, 0, 0, 0, -sin(q1), -cos(q1), 0],
                    [-cos(q1)*q3d, 0, 0, 0, 0, cos(q1), -sin(q1), 0],
                    [0, 4/5, 0, 0, 0, 0, 0, 6*q3d/5],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, -2*q3d, 0, 0]])
    B_sol = Matrix([])

    # Check that linearization is correct
    assert A.subs(upright_nominal) == A_sol
    assert B.subs(upright_nominal) == B_sol

    # Check eigenvalues at critical speed are all zero:
    assert A.subs(upright_nominal).subs(q3d, 1/sqrt(3)).eigenvals() == {0: 8}

def test_linearize_pendulum_minimal():
    q1 = dynamicsymbols('q1')                     # angle of pendulum
    u1 = dynamicsymbols('u1')                     # Angular velocity
    q1d = dynamicsymbols('q1', 1)                 # Angular velocity
    u1d = dynamicsymbols('u1', 1)                 # Angular acceleration
    L, m, t = symbols('L, m, t')
    g = 9.8

    # Compose world frame
    N = ReferenceFrame('N')
    pN = Point('N*')
    pN.set_vel(N, 0)

    # A.x is along the pendulum
    A = N.orientnew('A', 'axis', [q1, N.z])
    A.set_ang_vel(N, u1*N.z)

    # Locate point P relative to the origin N*
    P = pN.locatenew('P', L*A.x)
    P.v2pt_theory(pN, N, A)
    pP = Particle('pP', P, m)

    # Create Kinematic Differential Equations
    kde = Matrix([q1d - u1])

    # Input the force resultant at P
    R = m*g*N.x

    # Solve for eom with kanes method
    KM = KanesMethod(N, q_ind=[q1], u_ind=[u1], kd_eqs=kde)
    (fr, frstar) = KM.kanes_equations([(P, R)], [pP])

    # Linearize with Linearizer Class
    linearizer = KM.to_linearizer()
    A, B = linearizer.linearize(A_and_B=True)

    A_sol = Matrix([[0, 1], 
                    [-9.8*cos(q1)/L, 0]])
    B_sol = Matrix([])

    assert A == A_sol
    assert B == B_sol

def test_linearize_pendulum_nonminimal():
    # Create generalized coordinates and speeds for this non-minimal realization
    # q1, q2 = N.x and N.y coordinates of pendulum
    # u1, u2 = N.x and N.y velocities of pendulum
    q1, q2 = qs = dynamicsymbols('q1:3')
    q1d, q2d = qds = dynamicsymbols('q1:3', level=1)
    u1, u2 = us = dynamicsymbols('u1:3')
    u1d, u2d = uds = dynamicsymbols('u1:3', level=1)
    L, m, t = symbols('L, m, t')
    g = 9.8

    # Compose world frame
    N = ReferenceFrame('N')
    pN = Point('N*')
    pN.set_vel(N, 0)

    # A.x is along the pendulum
    theta1 = atan(q2/q1)
    A = N.orientnew('A', 'axis', [theta1, N.z])

    # Locate the pendulum mass
    P = pN.locatenew('P1', q1*N.x + q2*N.y)
    pP = Particle('pP', P, m)

    # Calculate the kinematic differential equations
    kde = Matrix([q1d - u1,
                  q2d - u2])
    dq_dict = solve(kde, [q1d, q2d])

    # Set velocity of point P
    P.set_vel(N, P.pos_from(pN).dt(N).subs(dq_dict))

    # Configuration constraint is length of pendulum
    f_c = Matrix([P.pos_from(pN).magnitude() - L])

    # Velocity constraint is that the velocity in the A.x direction is
    # always zero (the pendulum is never getting longer).
    f_v = Matrix([P.vel(N).express(A).dot(A.x)])
    f_v.simplify()

    # Acceleration constraints is the time derivative of the velocity constraint
    f_a = f_v.diff(t)
    f_a.simplify()

    # Input the force resultant at P
    R = m*g*N.x

    # Derive the equations of motion using the KanesMethod class.
    KM = KanesMethod(N, q_ind=[q1], u_ind=[u1], q_dependent=[q2],
            u_dependent=[u2], configuration_constraints=f_c,
            velocity_constraints=f_v, acceleration_constraints=f_a, kd_eqs=kde)
    (fr, frstar) = KM.kanes_equations([(P, R)], [pP])
    kdd = KM.kindiffdict()

    # Set the trim condition to be straight down, and non-moving
    q_trim = {q1: L, q2: 0}
    u_trim = {u1: 0, u2: 0}

    linearizer = KM.to_linearizer()
    A, B = linearizer.linearize(q_trim, u_trim, A_and_B=True)
    # TODO: This results in some nan and zoo (complex infinity) terms.
    # I don't think it should at this trim condition, but I could be wrong
