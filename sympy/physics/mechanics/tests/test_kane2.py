from sympy import cos, expand, Matrix, sin, symbols, tan, Poly, zeros, solve
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                     RigidBody, KanesMethod, inertia, Particle,
                                     dot, cross)

def test_aux_dep():
    # This test is about rolling disc dynamics, comparing the results found
    # with KanesMethod to those found when deriving the equations "manually"
    # with SymPy.
    # The terms Fr, Fr*, and Fr*_steady are all compared between the two
    # methods. Here, Fr*_steady refers to the generalized inertia forces for an
    # equilibrium configuration.
    # Note: comparing to the test of test_rolling_disc() in test_kane.py, this
    # test also tests auxiliary speeds and configuration and motion constraints
    #, seen in  the generalized dependent coordinates q[3], and depend speeds
    # u[3], u[4] and u[5].


    # First, mannual derivation of Fr, Fr_star, Fr_star_steady.

    # Symbols for time and constant parameters.
    # Symbols for contact forces: Fx, Fy, Fz.
    t, r, m, g, I, J = symbols('t r m g I J')
    Fx, Fy, Fz = symbols('Fx Fy Fz')

    # Configuration variables and their time derivatives:
    # q[0] -- yaw
    # q[1] -- lean
    # q[2] -- spin
    # q[3] -- dot(-r*B.z, A.z) -- distance from ground plane to disc center in
    #         A.z direction
    # Generalized speeds and their time derivatives:
    # u[0] -- disc angular velocity component, disc fixed x direction
    # u[1] -- disc angular velocity component, disc fixed y direction
    # u[2] -- disc angular velocity component, disc fixed z direction
    # u[3] -- disc velocity component, A.x direction
    # u[4] -- disc velocity component, A.y direction
    # u[5] -- disc velocity component, A.z direction
    # Auxiliary generalized speeds:
    # ua[0] -- contact point auxiliary generalized speed, A.x direction
    # ua[1] -- contact point auxiliary generalized speed, A.y direction
    # ua[2] -- contact point auxiliary generalized speed, A.z direction
    q = dynamicsymbols('q:4')
    qd = [qi.diff(t) for qi in q]
    u = dynamicsymbols('u:6')
    ud = [ui.diff(t) for ui in u]
    #ud_zero = {udi : 0 for udi in ud}
    ud_zero = dict(zip(ud, [0.]*len(ud)))
    ua = dynamicsymbols('ua:3')
    #ua_zero = {uai : 0 for uai in ua}
    ua_zero = dict(zip(ua, [0.]*len(ua)))

    # Reference frames:
    # Yaw intermediate frame: A.
    # Lean intermediate frame: B.
    # Disc fixed frame: C.
    N = ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', [q[0], N.z])
    B = A.orientnew('B', 'Axis', [q[1], A.x])
    C = B.orientnew('C', 'Axis', [q[2], B.y])

    # Angular velocity and angular acceleration of disc fixed frame
    # u[0], u[1] and u[2] are generalized independent speeds.
    C.set_ang_vel(N, u[0]*B.x + u[1]*B.y + u[2]*B.z)
    C.set_ang_acc(N, C.ang_vel_in(N).diff(t, B)
                   + cross(B.ang_vel_in(N), C.ang_vel_in(N)))

    # Velocity and acceleration of points:
    # Disc-ground contact point: P.
    # Center of disc: O, defined from point P with depend coordinate: q[3]
    # u[3], u[4] and u[5] are generalized dependent speeds.
    P = Point('P')
    P.set_vel(N, ua[0]*A.x + ua[1]*A.y + ua[2]*A.z)
    O = P.locatenew('O', q[3]*A.z + r*sin(q[1])*A.y)
    O.set_vel(N, u[3]*A.x + u[4]*A.y + u[5]*A.z)
    O.set_acc(N, O.vel(N).diff(t, A) + cross(A.ang_vel_in(N), O.vel(N)))

    # Kinematic differential equations:
    # Two equalities: one is w_c_n_qd = C.ang_vel_in(N) in three coordinates
    #                 directions of B, for qd0, qd1 and qd2.
    #                 the other is v_o_n_qd = O.vel(N) in A.z direction for qd3.
    # Then, solve for dq/dt's in terms of u's: qd_kd.
    w_c_n_qd = qd[0]*A.z + qd[1]*B.x + qd[2]*B.y
    v_o_n_qd = O.pos_from(P).diff(t, A) + cross(A.ang_vel_in(N), O.pos_from(P))
    kindiffs = Matrix([dot(w_c_n_qd - C.ang_vel_in(N), uv) for uv in B] +
                      [dot(v_o_n_qd - O.vel(N), A.z)])
    qd_kd = solve(kindiffs, qd)

    # Values of generalized speeds during a steady turn for later substitution
    # into the Fr_star_steady.
    steady_conditions = solve(kindiffs.subs({qd[1] : 0, qd[3] : 0}), u)
    steady_conditions.update({qd[1] : 0, qd[3] : 0})

    # Partial angular velocities and velocities.
    partial_w_C = [C.ang_vel_in(N).diff(ui, N) for ui in u + ua]
    partial_v_O = [O.vel(N).diff(ui, N) for ui in u + ua]
    partial_v_P = [P.vel(N).diff(ui, N) for ui in u + ua]

    # Configuration constraint: f_c, the projection of radius r in A.z direction
    #                                is q[3].
    # Velocity constraints: f_v, for u3, u4 and u5.
    # Acceleration constraints: f_a.
    f_c = Matrix([dot(-r*B.z, A.z) - q[3]])
    f_v = Matrix([dot(O.vel(N) - (P.vel(N) + cross(C.ang_vel_in(N),
        O.pos_from(P))), ai).expand() for ai in A])
    v_o_n = cross(C.ang_vel_in(N), O.pos_from(P))
    a_o_n = v_o_n.diff(t, A) + cross(A.ang_vel_in(N), v_o_n)
    f_a = Matrix([dot(O.acc(N) - a_o_n, ai) for ai in A])

    # Solve for constraint equations in the form of
    #                           u_dependent = A_rs * [u_i; u_aux].
    # First, obtain constraint coefficient matrix:  M_v * [u; ua] = 0;
    # Second, taking u[0], u[1], u[2] as independent,
    #         taking u[3], u[4], u[5] as dependent,
    #         rearranging the matrix of M_v to be A_rs for u_dependent.
    # Third, u_aux ==0 for u_dep, and resulting dictionary of u_dep_dict.
    M_v = zeros(3, 9)
    for i in range(3):
        for j, ui in enumerate(u + ua):
            M_v[i, j] = f_v[i].diff(ui)

    M_v_i = M_v[:, :3]
    M_v_d = M_v[:, 3:6]
    M_v_aux = M_v[:, 6:]
    M_v_i_aux = M_v_i.row_join(M_v_aux)
    A_rs = - M_v_d.inv() * M_v_i_aux

    u_dep = A_rs[:, :3] * Matrix(u[:3])
    u_dep_dict = dict(zip(u[3:], u_dep))
    #u_dep_dict = {udi : u_depi[0] for udi, u_depi in zip(u[3:], u_dep.tolist())}

    # Active forces: F_O acting on point O; F_P acting on point P.
    # Generalized active forces (unconstrained): Fr_u = F_point * pv_point.
    F_O = m*g*A.z
    F_P = Fx * A.x + Fy * A.y + Fz * A.z
    Fr_u = Matrix([dot(F_O, pv_o) + dot(F_P, pv_p) for pv_o, pv_p in
            zip(partial_v_O, partial_v_P)])

    # Inertia force: R_star_O.
    # Inertia of disc: I_C_O, where J is a inertia component about principal axis.
    # Inertia torque: T_star_C.
    # Generalized inertia forces (unconstrained): Fr_star_u.
    R_star_O = -m*O.acc(N)
    I_C_O = inertia(B, I, J, I)
    T_star_C = -(dot(I_C_O, C.ang_acc_in(N)) \
                 + cross(C.ang_vel_in(N), dot(I_C_O, C.ang_vel_in(N))))
    Fr_star_u = Matrix([dot(R_star_O, pv) + dot(T_star_C, pav) for pv, pav in
                        zip(partial_v_O, partial_w_C)])

    # Form nonholonomic Fr: Fr_c, and nonholonomic Fr_star: Fr_star_c.
    # Also, nonholonomic Fr_star in steady turning condition: Fr_star_steady.
    Fr_c = Fr_u[:3, :].col_join(Fr_u[6:, :]) + A_rs.T * Fr_u[3:6, :]
    Fr_star_c = Fr_star_u[:3, :].col_join(Fr_star_u[6:, :])\
                + A_rs.T * Fr_star_u[3:6, :]
    Fr_star_steady = Fr_star_c.subs(ud_zero).subs(u_dep_dict)\
            .subs(steady_conditions).subs({q[3]: -r*cos(q[1])}).expand()


    # Second, using KaneMethod in mechanics for fr, frstar and frstar_steady.

    # Rigid Bodies: disc, with inertia I_C_O.
    iner_tuple = (I_C_O, O)
    disc = RigidBody('disc', O, C, m, iner_tuple)
    bodyList = [disc]

    # Generalized forces: Gravity: F_o; Auxiliary forces: F_p.
    F_o = (O, F_O)
    F_p = (P, F_P)
    forceList = [F_o,  F_p]

    # KanesMethod.
    kane = KanesMethod(
        N, q_ind= q[:3], u_ind= u[:3], kd_eqs=kindiffs,
        q_dependent=q[3:], configuration_constraints = f_c,
        u_dependent=u[3:], velocity_constraints= f_v,
        u_auxiliary=ua
        )

    # fr, frstar, frstar_steady and kdd(kinematic differential equations).
    (fr, frstar)= kane.kanes_equations(forceList, bodyList)
    frstar_steady = frstar.subs(ud_zero).subs(u_dep_dict).subs(steady_conditions)\
                    .subs({q[3]: -r*cos(q[1])}).expand()
    kdd = kane.kindiffdict()


    # Test
    # First try Fr_c == fr;
    # Second try Fr_star_c == frstar;
    # Third try Fr_star_steady == frstar_steady.
    # Both signs are checked in case the equations were found with an inverse
    # sign.
    assert ((Matrix(Fr_c).expand() == fr.expand()) or
             (Matrix(Fr_c).expand() == (-fr).expand()))

    assert ((Matrix(Fr_star_c).expand() == frstar.expand()) or
             (Matrix(Fr_star_c).expand() == (-frstar).expand()))

    assert ((Matrix(Fr_star_steady).expand() == frstar_steady.expand()) or
             (Matrix(Fr_star_steady).expand() == (-frstar_steady).expand()))


def test_mat_inv_mul():
    # Just a quick test to check that KanesMethod._mat_inv_mul works as
    # intended. Uses SymPy generated primes as matrix entries, so each entry in
    # each matrix should be symbolic and unique, allowing proper comparison.
    # Checks _mat_inv_mul against Matrix.inv / Matrix.__mul__.
    from sympy import Matrix, prime
    from sympy.physics.mechanics import ReferenceFrame, KanesMethod

    # Just need to create an instance of KanesMethod to get to _mat_inv_mul
    mat_inv_mul = KanesMethod(ReferenceFrame('N'), [1], [1])._mat_inv_mul

    # going to form 3 matrices
    # 1 n x n
    # different n x n
    # 1 n x 2n
    n = 3
    m1 = Matrix(n, n, lambda i, j: prime(i * n + j + 2))
    m2 = Matrix(n, n, lambda i, j: prime(i * n + j + 5))
    m3 = Matrix(n, n, lambda i, j: prime(i + j * n + 2))

    assert mat_inv_mul(m1, m2) == m1.inv() * m2
    assert mat_inv_mul(m1, m3) == m1.inv() * m3
