import sympy as sm
import sympy.physics.mechanics as me

def test_removal_of_reaction_forces():
    """
    This function test whether the reaction forces are removed from the
    dynamic parts of (fr + frstar) and from the forcing vector.
    It is a simple pendulum where the foot is moving in the x direction.
    """
    t = me.dynamicsymbols._t
    N, A = me.ReferenceFrame('N'), me.ReferenceFrame('A')
    O, P, P1 = me.Point('O'), me.Point('P'), me.Point('P1')
    O.set_vel(N, 0)

    m, l, g = sm.symbols('m l g')
    q, x, y, u, ux, uy = me.dynamicsymbols('q x y u ux uy')
    x1, ux1 = me.dynamicsymbols('x1 u1')

    auxx, auxy, fx, fy, F = me.dynamicsymbols('auxx auxy fx fy F')

    A.orient_axis(N, q, N.z)
    A.set_ang_vel(N, u * N.z)
    P1.set_pos(O, x1 * N.x)
    P1.set_vel(N, ux*N.x + auxx*N.x + auxy*N.y)
    P.set_pos(P1, -l * A.y)
    P.v2pt_theory(P1, N, A)

    BODY = [me.Particle('body', P, m)]
    FL = [(P, -m * g * N.y), (P1, F*N.x + fx * N.x + fy * N.y)]

    kd = [u - q.diff(t), ux - x.diff(t), uy - y.diff(t), ux1 - x1.diff(t)]
    config_constr = [x - x1 - l * sm.sin(q), y - l * sm.cos(q)]
    speed_constr = [ux - ux1 - l * q.diff(t) * sm.cos(q),
        uy - l * q.diff(t) * sm.sin(q)]
    q_ind = [q, x1]
    q_dep = [x, y]
    u_ind = [u, ux1]
    u_dep = [ux, uy]
    aux = [auxx, auxy]

    reaction_forces = [fx, fy]

    KM = me.KanesMethod(
        N,
        q_ind=q_ind,
        q_dependent=q_dep,
        u_ind=u_ind,
        u_dependent=u_dep,
        u_auxiliary=aux,
        reaction_forces = reaction_forces,
        kd_eqs=kd,
        configuration_constraints=config_constr,
        velocity_constraints=speed_constr,
    )

    fr, frstar = KM.kanes_equations(BODY, FL)
    force = KM.forcing_full

    force_DS = me.find_dynamicsymbols(force)
    for i in reaction_forces + aux:
        assert i not in force_DS, f"{i} found in force_DS"

    eom_dynamic = sm.Matrix([(fr+frstar)[0], (fr+frstar)[1]])
    eom_DS = me.find_dynamicsymbols(eom_dynamic)
    for i in reaction_forces + aux:
        assert i not in eom_DS, f"{i} found in eom_DS"

    eom_dynamic = KM._form_eoms()
    eom_d = sm.Matrix([eom_dynamic[0], eom_dynamic[1]])
    eom_d_DS = me.find_dynamicsymbols(eom_d)
    for i in reaction_forces + aux:
        assert i not in eom_d_DS, f"{i} found in eom_d_DS"
