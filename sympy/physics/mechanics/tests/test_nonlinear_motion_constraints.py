from sympy import symbols, Matrix, solve
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
           Particle,  KanesMethod)
from sympy.testing import pytest
from sympy.solvers.solveset import NonlinearError

def test_nonlinear_motion_constraints():
    """
    This is the example in the paper "Forces associated with nonlinear,
    non-holonomic constrtaint equations", by Carlos M. Roithmaye and
    Dewey H. Hodges, Journal of Non-LinearApplied Mechanics 45 (2010),
    357 - 369
    It also verifies that a NonlinearError is raised when the non-linear
    velocity constraints are non-liner in udot.
    It also verifies that the reactions forces are same as the ones in the
    paper.

    """

    N = ReferenceFrame('N')
    O, P1, P2 = Point('O'), Point('P_1'), Point('P_2')

    q1, q2, q3, q4 = dynamicsymbols('q_1 q_2 q_3 q_4')
    u1, u2, u3, u4 = dynamicsymbols('u_1 u_2 u_3 u_4')

    aux1, aux2, aux3, aux4 = dynamicsymbols('aux1 aux2 aux3 aux4')
    f1, f2, f3, f4 = dynamicsymbols('f1 f2 f3 f4')
    t = dynamicsymbols._t
    O.set_vel(N, 0)

    sigma1, sigma2, sigma3, sigma4 = symbols('sigma_1 sigma_2 sigma_3 sigma_4')
    m1, m2 = symbols('m_1 m_2')

    P1.set_pos(O, q1 * N.x + q2 * N.y)
    P2.set_pos(O, q3 * N.x + q4 * N.y)
    P1.set_vel(N, u1*N.x + u2*N.y)
    P2.set_vel(N, u3*N.x + u4*N.y)

    P1a = Particle('P1a', P1, m1)
    P2a = Particle('P2a', P2, m2)
    bodies = [P1a, P2a]

    forces = [(P1, sigma1 * N.x + sigma2 * N.y), (P2, sigma3 * N.x + sigma4 * N.y)]


    nonlin_constr = Matrix([P1.vel(N).dot(P2.vel(N))])

    kd = Matrix([u1 - q1.diff(t), u2 - q2.diff(t), u3 - q3.diff(t), u4 - q4.diff(t)])

    q_ind = [q1, q2, q3, q4]
    u_ind = [u1, u2, u3]
    u_dep = [u4]
    aux = [aux1, aux2, aux3, aux4]

    KM = KanesMethod(
                N,
                q_ind,
                u_ind,
                u_dependent=u_dep,
                kd_eqs=kd,
                nonlinear_velocity_constraints=nonlin_constr,
    )

    fr, frstar = KM.kanes_equations(bodies, forces)

    #In the paper, uddot is solved for, using equation (33) und used when setting up
    #fr and frstar. As this is not done in Kane's method, I have to do it here.
    veldt = nonlin_constr.diff(t)
    repl = solve(veldt, [u4.diff(t)])
    frstar = (frstar.subs(repl))
    FR_sim = fr.col_join(frstar)

    # these are the equations (43), (44) and (45) from the paper
    fr1 = sigma1 - u3/u2 * sigma4
    fr2 = sigma2 - u4/u2 * sigma4
    fr3 = sigma3 - u1/u2 * sigma4

    # These are the original equations (47), (48) and (49) from the paper:
    #  frstar1 = -m1*u1.diff(t) - m2*u3/u2**2 * (u3*u1.diff(t) + u4*u2.diff(t) + u1*u3.diff(t))
    #  frstar2 = -m1*u2.diff(t) - m2*u4/u2**2 * (u3*u1.diff(t) + u4*u2.diff(t) + u1*u3.diff(t))
    #  frstar3 = -m2*u3.diff(t) - m2*u1/u2**2 * (u3*u1.diff(t) + u4*u2.diff(t) + u1*u3.diff(t))
    #
    # They have to be re-written as below to match the output of the Kane's method.

    frstar1 =-m1*u1.diff(t) + m2*(-u1*u3.diff(t)/u2 - u3*u1.diff(t)/u2 - u4*u2.diff(t)/u2)*u3/u2
    frstar2 = -m1*u2.diff(t) + m2*(-u1*u3.diff(t)/u2 - u3*u1.diff(t)/u2 - u4*u2.diff(t)/u2)*u4/u2
    frstar3 = m2*(-u1*u3.diff(t)/u2 - u3*u1.diff(t)/u2 - u4*u2.diff(t)/u2)*u1/u2 - m2*u3.diff(t)

    for i, j in enumerate((fr1, fr2, fr3, frstar1, frstar2, frstar3)):
        assert FR_sim[i] == j

    # This non-linear constraint is non-linear in udot
    nonlin_constr = Matrix([P1.vel(N).dot(P2.vel(N).diff(t, N))])

    with pytest.raises(NonlinearError):
        KanesMethod(N,
                    q_ind=q_ind,
                    u_ind=u_ind,
                    u_dependent=u_dep, kd_eqs=kd,
                    nonlinear_velocity_constraints=nonlin_constr
    )

    # Here the test for the reaction forces starts.
    P1.set_vel(N, u1*N.x + u2*N.y + aux1*N.x + aux2*N.y)
    P2.set_vel(N, u3*N.x + u4*N.y + aux3*N.x + aux4*N.y)
    forces = ([(P1, sigma1 * N.x + sigma2 * N.y), (P2, sigma3 * N.x + sigma4 * N.y)]
             + [(P1, f1*N.x + f2*N.y), (P2, f3*N.x + f4*N.y)])

    nonlin_constr = Matrix([P1.vel(N).dot(P2.vel(N))]).subs(dict.fromkeys([aux1, aux2, aux3, aux4], 0))

    KM = KanesMethod(
                N,
                q_ind,
                u_ind,
                u_dependent=u_dep,
                u_auxiliary=aux,
                kd_eqs=kd,
                nonlinear_velocity_constraints=nonlin_constr,
    )

    fr, frstar = KM.kanes_equations(bodies, forces)

    rhs = KM.rhs()
    rhs = rhs.subs(dict.fromkeys([f1, f2, f3, f4], 0))
    rhs_dict = {u1.diff(t): rhs[4], u2.diff(t): rhs[5], u3.diff(t): rhs[6], u4.diff(t): rhs[7]}
    reaction_forces = KM.auxiliary_eqs
    reaction_forces = (solve(reaction_forces, [f1, f2, f3, f4]))
    reaction_forces = {k: v.subs(rhs_dict) for k, v in reaction_forces.items()}


    # This is equation (34) of the a.m. article, just re-written so it matches
    # the output of the Kane's method
    lambda_34 = ((-m1*sigma3*u1 - m1*sigma4*u2 - m2*sigma1*u3 - m2*sigma2*u4) /
                 (m1*u1**2 + m1*u2**2 + m2*u3**2 + m2*u4**2))

    # This is the combination of eq(27) and eq(28) from the article
    C1x = lambda_34 * u3
    C1y = lambda_34 * u4
    C2x = lambda_34 * u1
    C2y = lambda_34 * u2
    results_article = [C1x, C1y, C2x, C2y]
    for i, j in enumerate(reaction_forces.keys()):
        test = (results_article[i] - reaction_forces[j]).simplify()
        assert test == 0
