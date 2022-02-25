from sympy.core.backend import (cos, sin, Matrix, symbols)
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                        KanesMethod, Particle, msubs)
from sympy.solvers.solvers import solve
from sympy.testing.pytest import raises

def test_replace_qdots_in_force():
    # Test PR 16700 "Replaces qdots with us in force-list in kanes.py"
    # The new functionality allows one to specify forces in qdots which will
    # automatically be replaced with u:s which are defined by the kde supplied
    # to KanesMethod. The test case is the double pendulum with interacting
    # forces in the example of chapter 4.7 "CONTRIBUTING INTERACTION FORCES"
    # in Ref. [1]. Reference list at end test function.

    q1, q2 = dynamicsymbols('q1, q2')
    qd1, qd2 = dynamicsymbols('q1, q2', level=1)
    u1, u2 = dynamicsymbols('u1, u2')

    l, m = symbols('l, m')

    N = ReferenceFrame('N') # Inertial frame
    A = N.orientnew('A', 'Axis', (q1, N.z)) # Rod A frame
    B = A.orientnew('B', 'Axis', (q2, N.z)) # Rod B frame

    O = Point('O') # Origo
    O.set_vel(N, 0)

    P = O.locatenew('P', ( l * A.x )) # Point @ end of rod A
    P.v2pt_theory(O, N, A)

    Q = P.locatenew('Q', ( l * B.x )) # Point @ end of rod B
    Q.v2pt_theory(P, N, B)

    Ap = Particle('Ap', P, m)
    Bp = Particle('Bp', Q, m)

    # Test will use different kdes
    kde1 = [u1 - qd1, u2 - qd2] # ref [1] formulation
    kde2 = [u1 - qd1, u2 - (qd1 + qd2)] # ref [2] formulation

    # The forces are specified below. sigma is the torsional spring stiffness
    # and delta is the viscous damping coefficient acting between the two
    # bodies. Here, we specify the viscous damper as function of qdots prior
    # forming the kde. In more complex systems it not might be obvious which
    # kde is most efficient, why it is convenient to specify viscous forces in
    # qdots independently of the kde.
    sig, delta = symbols('sigma, delta')

    # Will test both a linear and a cubic version of forces
    Ta = (sig * q2 + delta * qd2) * N.z
    forces_linear = [(A, Ta), (B, -Ta)]
    Ta_cubic = (sig * q2 + delta * qd2**3) * N.z
    forces_cubic  = [(A, Ta_cubic), (B, -Ta_cubic)]


    # Defining some helper functions
    def get_kde_sol(kdes):
        kde_sol = solve(kdes, [qd1, qd2])
        for _qd in [qd1, qd2]:
            kde_sol[_qd.diff()] = kde_sol[_qd].diff()
        return kde_sol

    def sub_qdots(forces, points, frames, kdes):
        # Sub qdots by u in forces, velocities, and accels
        # returns also a dictionary to restore velocities/accels
        # to using qdots in the expression
        kde_sol = get_kde_sol(kdes)
        
        original_qdots = {}
        for p in points:
            # before substituting, save original vel/acc as function of qdots
            original_qdots[p.set_vel] = (N, p.vel(N))
            original_qdots[p.set_acc] = (N, p.acc(N))
            
            p.set_vel(N, msubs(p.vel(N), kde_sol))
            p.set_acc(N, msubs(p.acc(N), kde_sol))
            
        for f in frames:
            original_qdots[f.set_ang_vel] = (N, f.ang_vel_in(N))
            original_qdots[f.set_ang_acc] = (N, f.ang_acc_in(N))
            
            f.set_ang_vel(N, msubs(f.ang_vel_in(N), kde_sol))
            f.set_ang_acc(N, msubs(f.ang_acc_in(N), kde_sol))
        
        forces_u = [(p, msubs(f, kde_sol)) for (p, f) in forces]
        return forces_u, original_qdots

    def restore_to_qdots(original_qdots):
        # restore velocities/accels to qdots
        for set_func, args in original_qdots.items():
            set_func(*args)
            
    for explicit_kinematics in [True, False]:
        original_qdots = {}
        for forces in [forces_linear, forces_cubic]:
            for kdi, kd_eqs in enumerate([kde1, kde2]):
                restore_to_qdots(original_qdots)

                KM = KanesMethod(N, [q1, q2], [u1, u2], kd_eqs=kd_eqs,
                    explicit_kinematics=explicit_kinematics)

                if explicit_kinematics:
                    fr, fstar = KM.kanes_equations([Ap, Bp], forces)
                else:
                    # We expect kanes_equations to complain about "qdot terms in velocities"
                    # when using implicit kinematics for this problem because
                    # it cant properly compute velocity partials
                    with raises(ValueError):
                        fr, fstar = KM.kanes_equations([Ap, Bp], forces)

                    # Try again but this time replace qdot terms with u 'manually'
                    # We do this to still ensure this problem gets tested although
                    # technically the main purpose of this test is to ensure
                    # that qdot terms in forces get replaced with u terms.        
                    forces_u, original_qdots = sub_qdots(forces, [P, Q], [A, B], kd_eqs)
                    fr, fstar = KM.kanes_equations([Ap, Bp], forces_u)

                if kdi == 0:
                    # Check fr with reference fr_expected from [1] with u:s instead of qdots.
                    if forces == forces_linear:
                        fr_expected = Matrix([ 0, -(sig*q2 + delta * u2) ])
                    else:
                        fr_expected = Matrix([ 0, -(sig*q2 + delta * u2**3) ])

                elif kdi == 1:            
                    if forces == forces_linear:
                        fr_expected = Matrix([sig * q2 + delta * (u2 - u1),
                                            - sig * q2 - delta * (u2 - u1)])
                    else:
                        fr_expected = Matrix([sig * q2 + delta * (u2 - u1)**3,
                                            - sig * q2 - delta * (u2 - u1)**3])

                    if forces == forces_linear:
                        # Only checking this for kdi==1 and Ta linear in qd1
                        # Mass and force matrix from p.6 in Ref. [2] with added forces from
                        # example of chapter 4.7 in [1] and without gravity.

                        forcing_matrix_expected = Matrix( [ [ m * l**2 * sin(q2) * u2**2 + sig * q2
                                                        + delta * (u2 - u1)],
                                                        [ m * l**2 * sin(q2) * -u1**2 - sig * q2
                                                        - delta * (u2 - u1)] ] )
                        mass_matrix_expected = Matrix( [ [ 2 * m * l**2, m * l**2 * cos(q2) ],
                                                    [ m * l**2 * cos(q2), m * l**2 ] ] )
                        assert (KM.mass_matrix.expand() == mass_matrix_expected.expand())
                        assert (KM.forcing.expand() == forcing_matrix_expected.expand())


                # This will be checked for both kd formulations and for linear and cubic forces
                assert fr.expand() == fr_expected.expand()

                # Specifying forces in u:s should stay the same
                # No need to test this when explicit_kinematics==False
                # as that's what's been tested already with sub_qdots above
                if explicit_kinematics:
                    kde_sol = get_kde_sol(kd_eqs)
                    forces_u = [(p, msubs(f, kde_sol)) for (p, f) in forces]
                    KM = KanesMethod(N, [q1, q2], [u1, u2], kd_eqs=kd_eqs, explicit_kinematics=explicit_kinematics)
                    fr, fstar = KM.kanes_equations([Ap, Bp], forces_u)

                    assert fr.expand() == fr_expected.expand()

        # References:
        # [1] T.R. Kane, D. a Levinson, Dynamics Theory and Applications, 2005.
        # [2] Arun K Banerjee, Flexible Multibody Dynamics:Efficient Formulations
        #     and Applications, John Wiley and Sons, Ltd, 2016.
        #     doi:http://dx.doi.org/10.1002/9781119015635.
