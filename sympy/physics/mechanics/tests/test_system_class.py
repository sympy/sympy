from sympy.core.backend import (Matrix, _simplify_matrix, zeros, sin, cos, eye,
                                symbols)
from sympy.simplify.simplify import simplify
from sympy.solvers.solvers import solve
from sympy.physics.mechanics import (
    dynamicsymbols, RigidBody, Particle, ReferenceFrame, System, PrismaticJoint,
    PinJoint, KanesMethod, Point, Torque, Force, LagrangesMethod)
from sympy.testing.pytest import raises

t = dynamicsymbols._t


def test_system_init():
    system = System()
    assert system.origin.name == 'inertial_origin'
    assert system.frame.name == 'inertial_frame'
    assert system.q_ind[:] == []
    assert system.q_dep[:] == []
    assert system.q[:] == []
    assert system.u_ind[:] == []
    assert system.u_dep[:] == []
    assert system.u[:] == []
    assert system.kdes[:] == []
    assert system.loads == tuple()
    assert system.bodies == tuple()
    assert system.joints == tuple()
    assert system.holonomic_constraints[:] == []
    assert system.nonholonomic_constraints[:] == []
    assert system.eom_method is None
    origin, frame = Point('O'), ReferenceFrame('N')
    system = System(origin, frame)
    assert system.origin == origin
    assert system.frame == frame
    rb = RigidBody('rb')
    system = system.from_newtonian(rb)
    assert system.bodies == (rb,)
    assert system.origin == rb.masscenter
    assert system.frame == rb.frame


def test_add_coordinates():
    q1, q2, q3, q4, q5, u1, u2, u3, u4, u5 = dynamicsymbols('q1:6 u1:6')
    system = System()
    # Test add_coordinates
    system.add_coordinates(q1, q2)
    assert system.q_ind == Matrix([q1, q2])
    assert system.q_dep[:] == []
    assert system.q == Matrix([q1, q2])
    system.add_coordinates(q3, independent=False)
    assert system.q_ind == Matrix([q1, q2])
    assert system.q_dep == Matrix([q3])
    assert system.q == Matrix([q1, q2, q3])
    system.add_coordinates(q4, q5, independent=False)
    assert system.q_ind == Matrix([q1, q2])
    assert system.q_dep == Matrix([q3, q4, q5])
    assert system.q == Matrix([q1, q2, q3, q4, q5])
    # Test add_speeds
    system.add_speeds(u1)
    assert system.u_ind == Matrix([u1])
    assert system.u_dep[:] == []
    assert system.u == Matrix([u1])
    system.add_speeds(u2, u3, independent=[True, False])
    assert system.u_ind == Matrix([u1, u2])
    assert system.u_dep == Matrix([u3])
    assert system.u == Matrix([u1, u2, u3])
    system.add_speeds(u4, u5, independent=[False, True])
    assert system.u_ind == Matrix([u1, u2, u5])
    assert system.u_dep == Matrix([u3, u4])
    assert system.u == Matrix([u1, u2, u5, u3, u4])
    system.q_ind = q1
    # Test setters
    assert system.q_ind == Matrix([q1])
    assert system.q_dep == Matrix([q3, q4, q5])
    assert system.q == Matrix([q1, q3, q4, q5])
    system.q_dep = [q2, q3]
    assert system.q_ind == Matrix([q1])
    assert system.q_dep == Matrix([q2, q3])
    assert system.q == Matrix([q1, q2, q3])
    system.u_ind = []
    assert system.u_ind[:] == []
    assert system.u_dep == Matrix([u3, u4])
    assert system.u == Matrix([u3, u4])
    system.u_dep = (u1, u2, u3)
    assert system.u_ind[:] == []
    assert system.u_dep == Matrix([u1, u2, u3])
    assert system.u == Matrix([u1, u2, u3])
    # Test if the col join works after resetting with a setter
    system.add_speeds(u4)
    assert system.u_ind == Matrix([u4])
    assert system.u_dep == Matrix([u1, u2, u3])
    assert system.u == Matrix([u4, u1, u2, u3])
    # Test incorrect data
    raises(ValueError, lambda: system.add_coordinates(q1))
    raises(ValueError, lambda: system.add_coordinates(q1, independent=False))
    raises(ValueError, lambda: system.add_speeds(q1, u5))
    raises(ValueError, lambda: system.add_speeds(u1))
    raises(ValueError, lambda: system.add_speeds(q1, independent=False))
    symb = symbols('symb')
    raises(ValueError, lambda: system.add_coordinates(symb))
    raises(ValueError, lambda: system.add_speeds(symb))
    with raises(ValueError):
        system.q_ind = [q1, q2]
    with raises(ValueError):
        system.q_dep = [u1]
    with raises(ValueError):
        system.u_ind = [u1, u2]
    with raises(ValueError):
        system.u_dep = [q3]
    # Nothing should have changed
    assert system.q_ind == Matrix([q1])
    assert system.q_dep == Matrix([q2, q3])
    assert system.q == Matrix([q1, q2, q3])
    assert system.u_ind == Matrix([u4])
    assert system.u_dep == Matrix([u1, u2, u3])
    assert system.u == Matrix([u4, u1, u2, u3])


def test_add_kdes():
    system = System()
    u1, u2, u3 = dynamicsymbols('u1:4')
    q1d, q2d, q3d = dynamicsymbols('q1:4', 1)
    # Test add_kdes
    system.add_kdes(q1d - u1)
    assert system.kdes == Matrix([q1d - u1])
    system.add_kdes(q2d + q3d - u2, q3d - q2d + u3)
    assert system.kdes == Matrix([q1d - u1, q2d + q3d - u2, q3d - q2d + u3])
    system.add_kdes(u1 - q1d)  # It should not recognize these
    assert system.kdes == Matrix([q1d - u1, q2d + q3d - u2, q3d - q2d + u3,
                                  u1 - q1d])
    # Test kdes setter
    system.kdes = q1d - u1
    assert system.kdes == Matrix([q1d - u1])
    system.kdes = [q1d - u1, u2 - q2d]
    assert system.kdes == Matrix([q1d - u1, u2 - q2d])
    system.kdes = []
    assert system.kdes[:] == []
    # Test if add_kdes still works after reset
    system.add_kdes(q1d - u1)
    assert system.kdes == Matrix([q1d - u1])
    # Test whether errors are raised
    raises(ValueError, lambda: system.add_kdes(q1d - u1))
    raises(TypeError, lambda: system.add_kdes([q3d - u3]))
    raises(ValueError, lambda: system.add_kdes(q1d - q1d))
    with raises(ValueError):
        system.kdes = [q1d - u1, u2 - q2d, q1d - u1]
    with raises(TypeError):
        system.kdes = 1
    assert system.kdes == Matrix([q1d - u1])


def test_add_constraints():
    q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1:4 u1:4')
    system = System()
    system.add_holonomic_constraints(q1 - q2, q3 - q1)
    # Test add methods
    assert system.holonomic_constraints == Matrix([q1 - q2, q3 - q1])
    assert system.nonholonomic_constraints[:] == []
    system.add_nonholonomic_constraints(u2 - u3)
    assert system.holonomic_constraints == Matrix([q1 - q2, q3 - q1])
    assert system.nonholonomic_constraints == Matrix([u2 - u3])
    # Test setters
    system.holonomic_constraints = q1 - q2
    assert system.holonomic_constraints == Matrix([q1 - q2])
    assert system.nonholonomic_constraints == Matrix([u2 - u3])
    system.nonholonomic_constraints = (u1 - 3, u2 - u1)
    assert system.holonomic_constraints == Matrix([q1 - q2])
    assert system.nonholonomic_constraints == Matrix([u1 - 3, u2 - u1])
    system.holonomic_constraints = []
    assert system.holonomic_constraints[:] == []
    # Test if add methods still works after reset
    system.add_holonomic_constraints(q3 - q1)
    assert system.holonomic_constraints == Matrix([q3 - q1])
    # Test whether errors are raised
    raises(ValueError, lambda: system.add_holonomic_constraints(q3 - q1))
    raises(TypeError, lambda: system.add_holonomic_constraints([q3 - q2]))
    raises(ValueError, lambda: system.add_holonomic_constraints(q1 - q1))
    with raises(ValueError):
        system.nonholonomic_constraints = [u2 - u1, u2 - u3, u2 - u1]
    with raises(TypeError):
        system.nonholonomic_constraints = 5 + 3
    assert system.holonomic_constraints == Matrix([q3 - q1])
    assert system.nonholonomic_constraints == Matrix([u1 - 3, u2 - u1])


def test_add_bodies():
    system = System()
    rb1, rb2 = RigidBody('rb1'), RigidBody('rb2')
    p1, p2 = Particle('p1'), Particle('p2')
    system.add_bodies(rb1, p1)
    assert system.bodies == (rb1, p1)
    system.add_bodies(p2)
    assert system.bodies == (rb1, p1, p2)
    system.bodies = []
    assert system.bodies == tuple()
    system.bodies = p2
    assert system.bodies == (p2,)
    symb = symbols('symb')
    raises(TypeError, lambda: system.add_bodies(symb))
    raises(ValueError, lambda: system.add_bodies(p2))
    with raises(TypeError):
        system.bodies = (rb1, rb2, p1, p2, symb)
    assert system.bodies == (p2,)


def test_add_loads():
    system = System()
    N, A = ReferenceFrame('N'), ReferenceFrame('A')
    rb1, rb2 = RigidBody('rb1', frame=N), RigidBody('rb2', frame=A)
    mc1, mc2 = Point('mc1'), Point('mc2')
    p1, p2 = Particle('p1', mc1), Particle('p2', mc2)
    system.add_loads(Torque(rb1, N.x), (mc1, A.x), Force(p1, A.x))
    assert system.loads == ((N, N.x), (mc1, A.x), (mc1, A.x))
    system.loads = [(A, A.x)]
    assert system.loads == ((A, A.x),)
    system.apply_force(rb1, N.x, rb2)
    assert system.loads == ((A, A.x), (rb1.masscenter, N.x),
                            (rb2.masscenter, -N.x))
    system.remove_load(rb1.masscenter)
    system.remove_load(rb2.masscenter)
    system.apply_torque(rb1, N.x, rb2)
    assert system.loads == ((A, A.x), (N, N.x), (A, -N.x))
    system.apply_force(p1, A.x)
    assert system.loads == ((A, A.x), (N, N.x), (A, -N.x), (mc1, A.x))
    raises(ValueError, lambda: system.add_loads((N, N.x, N.y)))
    with raises(TypeError):
        system.loads = (N, N.x)
    assert system.loads == ((A, A.x), (N, N.x), (A, -N.x), (mc1, A.x))
    system.clear_loads()
    assert system.loads == tuple()


def test_add_joints():
    q1, q2, q3, q4, u1, u2, u3 = dynamicsymbols('q1:5 u1:4')
    rb1, rb2, rb3, rb4, rb5 = symbols('rb1:6', cls=RigidBody)
    J1 = PinJoint('J1', rb1, rb2, q1, u1)
    J2 = PrismaticJoint('J2', rb2, rb3, q2, u2)
    J3 = PinJoint('J3', rb3, rb4, q3, u3)
    J_lag = PinJoint('J_lag', rb4, rb5, q4, q4.diff(t))
    system = System()
    system.add_joints(J1)
    assert system.joints == (J1,)
    assert system.bodies == (rb1, rb2)
    assert system.q_ind == Matrix([q1])
    assert system.u_ind == Matrix([u1])
    assert system.kdes == Matrix([u1 - q1.diff(t)])
    system.add_bodies(rb4)
    system.add_coordinates(q3)
    system.add_kdes(u3 - q3.diff(t))
    system.add_joints(J3)
    assert system.joints == (J1, J3)
    assert system.bodies == (rb1, rb2, rb4, rb3)
    assert system.q_ind == Matrix([q1, q3])
    assert system.u_ind == Matrix([u1, u3])
    assert system.kdes == Matrix([u1 - q1.diff(t), u3 - q3.diff(t)])
    system.add_joints(J2)
    assert system.joints == (J1, J3, J2)
    assert system.bodies == (rb1, rb2, rb4, rb3)
    assert system.q_ind == Matrix([q1, q3, q2])
    assert system.u_ind == Matrix([u1, u3, u2])
    assert system.kdes == Matrix([u1 - q1.diff(t), u3 - q3.diff(t),
                                  u2 - q2.diff(t)])
    system.add_joints(J_lag)
    assert system.joints == (J1, J3, J2, J_lag)
    assert system.bodies == (rb1, rb2, rb4, rb3, rb5)
    assert system.q_ind == Matrix([q1, q3, q2, q4])
    assert system.u_ind == Matrix([u1, u3, u2, q4.diff(t)])
    assert system.kdes == Matrix([u1 - q1.diff(t), u3 - q3.diff(t),
                                  u2 - q2.diff(t)])
    assert system.q_dep[:] == []
    assert system.u_dep[:] == []
    raises(ValueError, lambda: system.add_joints(J2))
    raises(TypeError, lambda: system.add_joints(rb1))


def test_get_body():
    system = System()
    rb1, rb2 = RigidBody('rb1'), RigidBody('rb2')
    p1, p2 = Particle('p1'), Particle('p2')
    system.add_bodies(rb1, rb2, p1, p2)
    assert system.get_body('p1') == p1
    assert system.get_body('rb2') == rb2
    assert system.get_body('rb3') is None


def test_get_joint():
    rb1, rb2, rb3, rb4 = symbols('rb1:5', cls=RigidBody)
    J1 = PinJoint('J1', rb1, rb2)
    J2 = PrismaticJoint('J2', rb2, rb3)
    J3 = PinJoint('J3', rb3, rb4)
    system = System()
    system.add_joints(J1, J2, J3)
    assert system.get_joint('J2') == J2
    assert system.get_joint('J1') == J1
    assert system.get_joint('J4') is None


def test_validate_system_kanes():
    # Create a correct full featured system
    q1, q2, q3, u1, u2, u3 = dynamicsymbols('q1:4 u1:4')
    q1d, q2d, q3d = dynamicsymbols('q1:4', 1)
    rb1, rb2, rb3, rb4 = symbols('rb1:5', cls=RigidBody)
    J1 = PinJoint('J1', rb1, rb2)
    J2 = PrismaticJoint('J2', rb2, rb3)
    J3 = PinJoint('J3', rb3, rb4)
    system = System()
    system.add_joints(J1, J2, J3)
    system.add_coordinates(q1, q2, q3, independent=[True, True, False])
    system.add_speeds(u1, u2, u3, independent=[True, False, False])
    system.add_kdes(q1d - u1, q2d + q3d - u2, q2d - q3d + u3)
    system.add_holonomic_constraints(q3 - q1 + q2)
    system.add_nonholonomic_constraints(u1 - q2d + u3)
    system.validate_system()
    q_ind, q_dep = system.q_ind[:], system.q_dep[:]
    u_ind, u_dep = system.u_ind[:], system.u_dep[:]
    kdes, hc = system.kdes[:], system.holonomic_constraints[:]
    nhc = system.nonholonomic_constraints[:]
    # Test Lagrange should fail due to the usage of generalized speeds
    raises(ValueError, lambda: system.validate_system(LagrangesMethod))
    # Test missing joint coordinate
    system.q_ind = q_ind[1:]
    system.u_ind = u_ind[:-1]
    system.kdes = kdes[:-1]
    raises(ValueError, lambda: system.validate_system())
    # Test missing joint speed
    system.q_ind = q_ind[:-1]
    system.u_ind = u_ind[1:]
    raises(ValueError, lambda: system.validate_system())
    # Test missing joint kdes
    system.u_ind = u_ind[:-1]
    system.validate_system()
    system.kdes = kdes[1:]
    raises(ValueError, lambda: system.validate_system())
    # Test missing holonomic constraint
    system.q_ind = q_ind
    system.u_ind = u_ind
    system.kdes = kdes
    system.validate_system()
    system.holonomic_constraints = []
    system.nonholonomic_constraints = nhc + [u3 - u1 + u2]
    raises(ValueError, lambda: system.validate_system())
    system.q_dep = []
    system.q_ind = q_ind + [q3]
    system.validate_system()
    # Test missing nonholonomic constraint
    system.q_ind = q_ind
    system.q_dep = q_dep
    system.holonomic_constraints = hc
    system.nonholonomic_constraints = []
    raises(ValueError, lambda: system.validate_system())
    system.u_dep = u2
    system.u_ind = u_ind + [u3]
    system.validate_system()
    # Test more speeds than coordinates
    system.q_ind = q_ind[:-1]
    system.u_ind = u_ind
    system.u_dep = u_dep
    system.nonholonomic_constraints = nhc
    system.validate_system()
    # Test more coordinates than speeds
    system.q_ind = q_ind
    system.u_ind = u_ind[:-1]
    system.kdes = kdes[:-1]
    raises(ValueError, lambda: system.validate_system())
    # Test wrong number of kdes
    system.u_ind = u_ind
    raises(ValueError, lambda: system.validate_system())
    system.kdes = kdes + [u3 + u2 - q3d]
    raises(ValueError, lambda: system.validate_system())
    # Test if check_duplicates does not result in errors
    system.kdes = kdes
    system.validate_system(check_duplicates=True)


def test_validate_system_lagrange():
    # Create a correct full featured system
    qj1, qj2, qj3, q1, q2, q3, u = dynamicsymbols('qj1:4 q1:4 u')
    qj1d, qj2d, qj3d, q1d, q2d, q3d = dynamicsymbols('qj1:4 q1:4', 1)
    rb1, rb2, rb3, rb4 = symbols('rb1:5', cls=RigidBody)
    J1 = PinJoint('J1', rb1, rb2, qj1, qj1d)
    J2 = PrismaticJoint('J2', rb2, rb3, qj2, qj2d)
    J3 = PinJoint('J3', rb3, rb4, qj3, qj3d)
    system = System()
    system.add_joints(J1, J2, J3)
    system.add_coordinates(q1, q2, q3, independent=[True, True, False])
    system.add_holonomic_constraints(q3 - q1 + q2)
    system.add_nonholonomic_constraints(q1d - q2d + q3d)
    system.validate_system(LagrangesMethod)
    system.u_ind = []  # Should not matter
    system.validate_system(LagrangesMethod)
    q_ind, q_dep = system.q_ind[:], system.q_dep[:]
    hc = system.holonomic_constraints[:]
    nhc = system.nonholonomic_constraints[:]
    # KanesMethod should throw errors
    raises(ValueError, lambda: system.validate_system())
    # Test missing joint coordinate
    system.q_ind = q_ind[1:]
    raises(ValueError, lambda: system.validate_system(LagrangesMethod))
    # Test missing holonomic constraint
    system.q_ind = q_ind
    system.holonomic_constraints = []
    system.nonholonomic_constraints = nhc + [q3d - q1d + q2d]
    raises(ValueError, lambda: system.validate_system(LagrangesMethod))
    system.q_dep = []
    system.q_ind = q_ind + [q3]
    system.validate_system(LagrangesMethod)
    # Test adding generalized speeds
    system.q_ind = q_ind
    system.q_dep = q_dep
    system.holonomic_constraints = hc
    system.nonholonomic_constraints = []
    system.validate_system(LagrangesMethod)
    system.u_ind = u
    raises(ValueError, lambda: system.validate_system(LagrangesMethod))
    system.u_ind = []
    system.u_dep = u
    raises(ValueError, lambda: system.validate_system(LagrangesMethod))


def test_cart_pendulum_kanes():
    # This example is the same as in the top documentation of System
    g, l, mc, mp = symbols('g l mc mp')
    F, qp, qc, up, uc = dynamicsymbols('F qp qc up uc')
    rail = RigidBody('rail')
    cart = RigidBody('cart', mass=mc)
    bob = Particle('bob', mass=mp)
    bob_frame = ReferenceFrame('bob_frame')
    system = System.from_newtonian(rail)
    assert system.bodies == (rail,)
    assert system.frame == rail.frame
    assert system.origin == rail.masscenter
    slider = PrismaticJoint('slider', rail, cart, qc, uc, joint_axis=rail.x)
    pin = PinJoint('pin', cart, bob, qp, up, joint_axis=cart.z,
                   child_interframe=bob_frame, child_point=l * bob_frame.y)
    system.add_joints(slider, pin)
    assert system.joints == (slider, pin)
    assert system.get_joint('slider') == slider
    assert system.get_body('bob') == bob
    system.apply_gravity(-g * system.y)
    system.apply_force(cart, F * rail.x)
    system.validate_system()
    system.form_eoms()
    assert isinstance(system.eom_method, KanesMethod)
    assert (_simplify_matrix(system.mass_matrix - Matrix(
        [[mp + mc, mp * l * cos(qp)], [mp * l * cos(qp), mp * l ** 2]]))
            == zeros(2, 2))
    assert (_simplify_matrix(system.forcing - Matrix([
        [mp * l * up ** 2 * sin(qp) + F], [-mp * g * l * sin(qp)]]))
            == zeros(2, 1))

    system.add_holonomic_constraints(
        bob.masscenter.pos_from(rail.masscenter).dot(system.x))
    assert system.eom_method is None
    system.q_ind, system.q_dep = qp, qc
    system.u_ind, system.u_dep = up, uc
    system.validate_system()

    # Computed solution based on manually solving the constraints
    subs = {qc: -l * sin(qp),
            uc: -l * cos(qp) * up,
            uc.diff(t): l * (up ** 2 * sin(qp) - up.diff(t) * cos(qp))}
    upd_expected = ((-g * mp * sin(qp) + l * mc * sin(2 * qp) * up ** 2 / 2 -
                     l * mp * sin(2 * qp) * up ** 2 / 2 - F * cos(qp)) /
                    (l * (mc * cos(qp) ** 2 + mp * sin(qp) ** 2)))
    upd_sol = solve(system.form_eoms().xreplace(subs), up.diff(t))[up.diff(t)]
    assert simplify(upd_sol - upd_expected) == 0
    assert isinstance(system.eom_method, KanesMethod)

    # Test other output
    Mk = -Matrix([[0, 1], [1, 0]])
    gk = -Matrix([uc, up])
    Md = Matrix([[-l ** 2 * mp * cos(qp) ** 2 + l ** 2 * mp,
                  l * mp * cos(qp) - l * (mc + mp) * cos(qp)],
                 [l * cos(qp), 1]])
    gd = Matrix([[-g * l * mp * sin(qp) - l ** 2 * mp * up ** 2 * sin(qp) * cos(
        qp) - l * F * cos(qp)], [l * up ** 2 * sin(qp)]])
    Mm = (Mk.row_join(zeros(2, 2))).col_join(zeros(2, 2).row_join(Md))
    gm = gk.col_join(gd)
    assert _simplify_matrix(system.mass_matrix - Md) == zeros(2, 2)
    assert _simplify_matrix(system.forcing - gd) == zeros(2, 1)
    assert _simplify_matrix(system.mass_matrix_full - Mm) == zeros(4, 4)
    assert _simplify_matrix(system.forcing_full - gm) == zeros(4, 1)


def test_cart_pendulum_lagrange():
    # Lagrange version of test_cart_pendulus_kanes
    g, l, mc, mp = symbols('g l mc mp')
    F, qp, qc = dynamicsymbols('F qp qc')
    qpd, qcd = dynamicsymbols('qp qc', 1)
    rail = RigidBody('rail')
    cart = RigidBody('cart', mass=mc)
    bob = Particle('bob', mass=mp)
    bob_frame = ReferenceFrame('bob_frame')
    system = System.from_newtonian(rail)
    assert system.bodies == (rail,)
    assert system.frame == rail.frame
    assert system.origin == rail.masscenter
    slider = PrismaticJoint('slider', rail, cart, qc, qcd, joint_axis=rail.x)
    pin = PinJoint('pin', cart, bob, qp, qpd, joint_axis=cart.z,
                   child_interframe=bob_frame, child_point=l * bob_frame.y)
    system.add_joints(slider, pin)
    assert system.joints == (slider, pin)
    assert system.get_joint('slider') == slider
    assert system.get_body('bob') == bob
    for body in system.bodies:
        body.potential_energy = body.mass * g * body.masscenter.pos_from(
            system.origin).dot(system.y)
    system.apply_force(cart, F * rail.x)
    system.validate_system(LagrangesMethod)
    system.form_eoms(LagrangesMethod)
    assert (_simplify_matrix(system.mass_matrix - Matrix(
        [[mp + mc, mp * l * cos(qp)], [mp * l * cos(qp), mp * l ** 2]]))
            == zeros(2, 2))
    assert (_simplify_matrix(system.forcing - Matrix([
        [mp * l * qpd ** 2 * sin(qp) + F], [-mp * g * l * sin(qp)]]))
            == zeros(2, 1))

    system.add_holonomic_constraints(
        bob.masscenter.pos_from(rail.masscenter).dot(system.x))
    assert system.eom_method is None
    system.q_ind, system.q_dep = qp, qc

    # Computed solution based on manually solving the constraints
    subs = {qc: -l * sin(qp),
            qcd: -l * cos(qp) * qpd,
            qcd.diff(t): l * (qpd ** 2 * sin(qp) - qpd.diff(t) * cos(qp))}
    qpdd_expected = ((-g * mp * sin(qp) + l * mc * sin(2 * qp) * qpd ** 2 / 2 -
                      l * mp * sin(2 * qp) * qpd ** 2 / 2 - F * cos(qp)) /
                     (l * (mc * cos(qp) ** 2 + mp * sin(qp) ** 2)))
    eoms = system.form_eoms(LagrangesMethod)
    lam1 = system.eom_method.lam_vec[0]
    lam1_sol = system.eom_method.solve_multipliers()[lam1]
    qpdd_sol = solve(eoms[0].xreplace({lam1: lam1_sol}).xreplace(subs),
                     qpd.diff(t))[0]
    assert simplify(qpdd_sol - qpdd_expected) == 0
    assert isinstance(system.eom_method, LagrangesMethod)

    # Test other output
    Md = Matrix([[l ** 2 * mp, l * mp * cos(qp), -l * cos(qp)],
                 [l * mp * cos(qp), mc + mp, -1]])
    gd = Matrix([[-g * l * mp * sin(qp)], [l * mp * sin(qp) * qpd ** 2 + F]])
    Mm = (eye(2).row_join(zeros(2, 3))).col_join(zeros(3, 2).row_join(
        Md.col_join(Matrix([l * cos(qp), 1, 0]).T)))
    gm = Matrix([qpd, qcd] + gd[:] + [l * sin(qp) * qpd ** 2])
    assert _simplify_matrix(system.mass_matrix - Md) == zeros(2, 3)
    assert _simplify_matrix(system.forcing - gd) == zeros(2, 1)
    assert _simplify_matrix(system.mass_matrix_full - Mm) == zeros(5, 5)
    assert _simplify_matrix(system.forcing_full - gm) == zeros(5, 1)
