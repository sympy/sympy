import sympy as sm
from sympy.external import import_module
import sympy.physics.mechanics as me
import sympy.physics._biomechanics as bm
from sympy.testing.pytest import skip
numpy = import_module('numpy')
scipy = import_module('scipy', import_kwargs={'fromlist': ['optimize']})


def test_arm_lever_tutorial():
    if not numpy:
        skip('NumPy not installed.')
    if not scipy:
        skip('SciPy not installed.')

    q1, q2, q3, q4 = me.dynamicsymbols('q1, q2, q3, q4')
    u1, u2, u3, u4 = me.dynamicsymbols('u1, u2, u3, u4')

    dx, dy, dz, lA, lC, lD = sm.symbols('dx, dy, dz, lA, lC, lD', real=True,
                                        nonnegative=True)
    mA, mC, mD = sm.symbols('mA, mC, mD', real=True, positive=True)
    g, k, c, r = sm.symbols('g, k, c, r', real=True, positive=True)

    N, A, B, C, D = sm.symbols('N, A, B, C, D', cls=me.ReferenceFrame)
    O, P1, P2, P3, P4 = sm.symbols('O, P1, P2, P3, P4 ', cls=me.Point)
    Ao, Co, Cm, Dm, Do = sm.symbols('Ao, Co, Cm, Dm, Do', cls=me.Point)

    A.orient_axis(N, q1, N.z)
    B.orient_axis(N, q2, N.y)
    C.orient_axis(B, q3, B.z)
    D.orient_axis(C, q4, C.y)
    A.set_ang_vel(N, u1*N.z)
    B.set_ang_vel(N, u2*N.y)
    C.set_ang_vel(B, u3*B.z)
    D.set_ang_vel(C, u4*C.y)

    Ao.set_pos(O, dx*N.x)
    P1.set_pos(Ao, lA*A.y)
    P2.set_pos(O, dy*N.y + dz*N.z)
    Co.set_pos(P2, lC/2*C.z)
    Cm.set_pos(P2, 1*lC/3*C.z)
    P3.set_pos(P2, lC*C.z)
    Dm.set_pos(P3, 1*lD/3*D.z)
    Do.set_pos(P3, lD/2*D.z)
    P4.set_pos(P3, lD*D.z)

    O.set_vel(N, 0)
    Ao.set_vel(N, 0)
    P1.v2pt_theory(Ao, N, A)
    P2.set_vel(N, 0)
    Co.v2pt_theory(P2, N, C)
    Cm.v2pt_theory(P2, N, C)
    P3.v2pt_theory(P2, N, C)
    Dm.v2pt_theory(P3, N, D)
    Do.v2pt_theory(P3, N, D)
    P4.v2pt_theory(P3, N, D)

    holonomic = (P4.pos_from(O) - P1.pos_from(O)).to_matrix(N)

    IA = me.Inertia(me.inertia(A, mA/12*lA**2, mA/2*lA**2, mA/12*lA**2), Ao)
    IC = me.Inertia(me.inertia(C, mC/12*lC**2, mC/12*lC**2, mC/2*lC**2), Co)
    ID = me.Inertia(me.inertia(D, mD/12*lD**2, mD/12*lD**2, mD/2*lD**2), Do)

    lever = me.RigidBody('lever', masscenter=Ao, frame=A, mass=mA, inertia=IA)
    u_arm = me.RigidBody('upper arm', masscenter=Co, frame=C, mass=mC,
                         inertia=IC)
    l_arm = me.RigidBody('lower arm', masscenter=Do, frame=D, mass=mD,
                         inertia=ID)

    lever_resistance = me.Torque(A, (-k*q1 - c*u1)*N.z)

    gravC = me.Force(u_arm, mC*g*N.z)
    gravD = me.Force(l_arm, mD*g*N.z)

    bicep_pathway = me.LinearPathway(Cm, Dm)

    bicep_activation = \
        bm.FirstOrderActivationDeGroote2016.with_default_constants('bicep')

    bicep = bm.MusculotendonDeGroote2016('bicep', bicep_pathway,
                                         activation_dynamics=bicep_activation)

    class ExtensorPathway(me.PathwayBase):

        def __init__(self, origin, insertion, axis_point, axis, parent_axis,
                     child_axis, radius, coordinate):
            """A custom pathway that wraps a cicular arc around a pin joint.

            This is intended to be used for extensor muscles. For example, a
            tricep wrapping around the elbow joint to extend the upper arm at
            the elbow.

            Parameters
            ==========
            origin : Point
                Muscle origin point fixed on the parent body (A).
            insertion : Point
                Muscle insertion point fixed on the child body (B).
            axis_point : Point
                Pin joint location fixed in both the parent and child.
            axis : Vector
                Pin joint rotation axis.
            parent_axis : Vector
                Axis fixed in the parent frame (A) that is directed from the
                pin joint point to the muscle origin point.
            child_axis : Vector
                Axis fixed in the child frame (B) that is directed from the pin
                joint point to the muscle insertion point.
            radius : sympyfiable
                Radius of the arc that the muscle wraps around.
            coordinate : sympfiable function of time
                Joint angle, zero when parent and child frames align. Positive
                rotation about the pin joint axis, B with respect to A.

            Notes
            =====

            Only valid for coordinate >= 0.

            """
            super().__init__(origin, insertion)

            self.origin = origin
            self.insertion = insertion
            self.axis_point = axis_point
            self.axis = axis.normalize()
            self.parent_axis = parent_axis.normalize()
            self.child_axis = child_axis.normalize()
            self.radius = radius
            self.coordinate = coordinate

            self.origin_distance = axis_point.pos_from(origin).magnitude()
            self.insertion_distance = axis_point.pos_from(insertion).magnitude()
            self.origin_angle = sm.asin(self.radius/self.origin_distance)
            self.insertion_angle = sm.asin(self.radius/self.insertion_distance)

        @property
        def length(self):
            """Length of the pathway.

            Length of two fixed length line segments and a changing arc length
            of a circle.

            """

            angle = self.origin_angle + self.coordinate + self.insertion_angle
            arc_length = self.radius*angle

            origin_segment_length = (self.origin_distance *
                                     sm.cos(self.origin_angle))
            insertion_segment_length = (self.insertion_distance *
                                        sm.cos(self.insertion_angle))

            return origin_segment_length + arc_length + insertion_segment_length

        @property
        def extension_velocity(self):
            """Extension velocity of the pathway.

            Arc length of circle is the only thing that changes when the elbow
            flexes and extends.

            """
            return self.radius*self.coordinate.diff(me.dynamicsymbols._t)

        def compute_loads(self, force_magnitude):
            """Loads in the correct format to be supplied to `KanesMethod`.

            Forces applied to origin, insertion, and P from the muscle wrapped
            over circular arc of radius r.

            """

            parent_tangency_point = me.Point('Aw')  # fixed in parent
            child_tangency_point = me.Point('Bw')  # fixed in child

            parent_tangency_point.set_pos(
                self.axis_point,
                -(self.radius*sm.cos(self.origin_angle) *
                  self.parent_axis.cross(self.axis))
                + self.radius*sm.sin(self.origin_angle)*self.parent_axis,
            )
            child_tangency_point.set_pos(
                self.axis_point,
                (self.radius*sm.cos(self.insertion_angle) *
                 self.child_axis.cross(self.axis))
                + self.radius*sm.sin(self.insertion_angle)*self.child_axis),

            parent_force_direction_vector = self.origin.pos_from(parent_tangency_point)
            child_force_direction_vector = self.insertion.pos_from(child_tangency_point)
            force_on_parent = force_magnitude*parent_force_direction_vector.normalize()
            force_on_child = force_magnitude*child_force_direction_vector.normalize()
            loads = [
                me.Force(self.origin, force_on_parent),
                me.Force(self.axis_point, -(force_on_parent + force_on_child)),
                me.Force(self.insertion, force_on_child),
            ]
            return loads

    tricep_pathway = ExtensorPathway(Cm, Dm, P3, B.y, -C.z, D.z, r, q4)
    tricep_activation = \
        bm.FirstOrderActivationDeGroote2016.with_default_constants('tricep')
    tricep = bm.MusculotendonDeGroote2016('tricep', tricep_pathway,
                                          activation_dynamics=tricep_activation)

    loads = (
        bicep.to_loads() +
        tricep.to_loads() +
        [lever_resistance, gravC, gravD]
    )

    kane = me.KanesMethod(
        N,
        (q1,),
        (u1,),
        kd_eqs=(
            u1 - q1.diff(),
            u2 - q2.diff(),
            u3 - q3.diff(),
            u4 - q4.diff(),
        ),
        q_dependent=(q2, q3, q4),
        configuration_constraints=holonomic,
        velocity_constraints=holonomic.diff(me.dynamicsymbols._t),
        u_dependent=(u2, u3, u4),
    )
    Fr, Frs = kane.kanes_equations((lever, u_arm, l_arm), loads)

    dadt = bicep.rhs().col_join(tricep.rhs())

    q, u = kane.q, kane.u
    a = bicep.x.col_join(tricep.x)
    x = q.col_join(u).col_join(a)
    e = bicep.r.col_join(tricep.r)

    p = sm.Matrix([
        dx,
        dy,
        dz,
        lA,
        lC,
        lD,
        mA,
        mC,
        mD,
        g,
        k,
        c,
        r,
        bicep._F_M_max,
        bicep._l_M_opt,
        bicep._l_T_slack,
        bicep._v_M_max,
        bicep._alpha_opt,
        bicep._beta,
        tricep._F_M_max,
        tricep._l_M_opt,
        tricep._l_T_slack,
        tricep._v_M_max,
        tricep._alpha_opt,
        tricep._beta,
    ])
    p

    eval_diffeq = sm.lambdify((q, u, a, e, p),
                              (kane.mass_matrix, kane.forcing, dadt), cse=True)
    eval_holonomic = sm.lambdify((q, p), holonomic, cse=True)

    p_vals = numpy.array([
        0.31,  # dx [m]
        0.15,  # dy [m]
        -0.31,  # dz [m]
        0.2,   # lA [m]
        0.3,  # lC [m]
        0.3,  # lD [m]
        0.5,  # mA [kg]
        2.3,  # mC [kg]
        1.7,  # mD [kg]
        9.81,  # g [m/s/s]
        5.0,  # k [Nm/rad]
        0.5,  # c [Nms/rad]
        0.03,  # r [m]
        500.0,
        0.6*0.3,
        0.55*0.3,
        10.0,
        0.0,
        0.1,
        500.0,
        0.6*0.3,
        0.65*0.3,
        10.0,
        0.0,
        0.1,
    ])

    q_vals = numpy.array([
        numpy.deg2rad(5.0),  # q1 [rad]
        numpy.deg2rad(-10.0),  # q2 [rad]
        numpy.deg2rad(0.0),  # q3 [rad]
        numpy.deg2rad(75.0),  # q4 [rad]
    ])

    def eval_holo_fsolve(x):
        q1 = q_vals[0]  # specified
        q2, q3, q4 = x
        return eval_holonomic((q1, q2, q3, q4), p_vals).squeeze()

    q_vals[1:] = scipy.optimize.fsolve(eval_holo_fsolve, q_vals[1:])

    u_vals = numpy.array([
        0.0,  # u1, [rad/s]
        0.0,  # u2, [rad/s]
        0.0,  # u3, [rad/s]
        0.0,  # u4, [rad/s]
    ])

    a_vals = numpy.array([
        0.0,  # a_bicep, nondimensional
        0.0,  # a_tricep, nondimensional
    ])

    e_vals = numpy.array([
        0.0,
        0.0,
    ])

    def eval_rhs(t, x, p):

        q = x[0:4]
        u = x[4:8]
        a = x[8:10]

        if t < 0.5 or t > 1.5:
            e = numpy.array([0.0, 0.0])
        else:
            e = numpy.array([0.8, 0.0])

        qd = u
        m, f, ad = eval_diffeq(q, u, a, e, p)
        ud = numpy.linalg.solve(m, f).squeeze()

        return numpy.hstack((qd, ud, ad.squeeze()))

    x0 = numpy.hstack((q_vals, u_vals, a_vals))

    expected_rhs = numpy.array([0.0, 0.0, 0.0, 0.0, -8.86014888, 5.96818164,
                                0.54422603, -5.82455426, 0.0, 0.0])
    numpy.testing.assert_allclose(expected_rhs, eval_rhs(0.0, x0, p_vals))
