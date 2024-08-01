def test_linearity_of_motion_constraints():
    # Test that an error is raised if nonlinear velocity constraints are
    # supplied to the KanesMethod.
    # It is a simple pendulum.
    t =dynamicsymbols._t
    N, A = ReferenceFrame('N'), ReferenceFrame('A')
    O, P = Point('O'), Point('P')
    O.set_vel(N, 0)

    m, l, g = symbols('m l g')
    q, x, y, u, ux, uy = dynamicsymbols('q x y u ux uy')

    A.orient_axis(N, q, N.z)
    A.set_ang_vel(N, u * N.z)
    P.set_pos(O, -l * A.y)
    P.v2pt_theory(O, N, A)

    kd = [u - q.diff(t), ux - x.diff(t), uy - y.diff(t)]
    config_constr = [x - l * sin(q), y - l * cos(q)]

    q_ind = [q]
    q_dep = [x, y]
    u_ind = [u]
    u_dep = [ux, uy]

    # Make sure an error is raised if nonlinear velocity constraints are
    # supplied.
    speed_constr = [ux**2 - l * q.diff(t) * cos(q),sin(uy) +
        l * q.diff(t)**2 * sin(q)]

    KM = KanesMethod(
        N,
        q_ind=q_ind,
        q_dependent=q_dep,
        u_ind=u_ind,
        u_dependent=u_dep,
        kd_eqs=kd,
        configuration_constraints=config_constr,
        velocity_constraints=speed_constr,
    )

    raises(ValueError, lambda: KM)