import pytest

from sympy.core.backend import ImmutableMatrix, symbols
from sympy.physics.mechanics import (
    System, ReferenceFrame, Point, RigidBody, Particle, dynamicsymbols,
    PinJoint, PrismaticJoint, Force, Torque)

t = dynamicsymbols._t  # type: ignore
q = dynamicsymbols('q:6')  # type: ignore
qd = dynamicsymbols('q:6', 1)  # type: ignore
u = dynamicsymbols('u:6')  # type: ignore


class TestSystem:
    @pytest.fixture()
    def _empty_system_setup(self):
        self.system = System(Point('origin'), ReferenceFrame('frame'))

    def _empty_system_check(self, exclude=()):
        matrices = ('q_ind', 'q_dep', 'q', 'u_ind', 'u_dep', 'u', 'kdes',
                    'holonomic_constraints', 'nonholonomic_constraints')
        tuples = ('loads', 'bodies', 'joints')
        for attr in matrices:
            if attr not in exclude:
                assert getattr(self.system, attr)[:] == []
        for attr in tuples:
            if attr not in exclude:
                assert getattr(self.system, attr) == ()
        if 'eom_method' not in exclude:
            assert self.system.eom_method is None

    def test_empty_system(self, _empty_system_setup):
        self._empty_system_check()

    @pytest.fixture()
    def _filled_system_setup(self, _empty_system_setup):
        self.bodies = symbols('rb1:5', cls=RigidBody)
        self.joints = (
            PinJoint('J1', self.bodies[0], self.bodies[1], q[0], u[0]),
            PrismaticJoint('J2', self.bodies[1], self.bodies[2], q[1], u[1]),
            PinJoint('J3', self.bodies[2], self.bodies[3], q[2], u[2])
        )
        self.system.add_joints(*self.joints)
        self.system.add_coordinates(q[3], independent=[False])
        self.system.add_speeds(u[3], independent=False)
        self.system.add_kdes(u[3] - qd[3])
        self.system.add_holonomic_constraints(q[2] - q[0] + q[1])
        self.system.add_nonholonomic_constraints(u[3] - qd[1] + u[2])
        self.system.u_ind = u[:2]
        self.system.u_dep = u[2:4]

    def _filled_system_check(self, exclude=()):
        assert 'q_ind' in exclude or self.system.q_ind[:] == q[:3]
        assert 'q_dep' in exclude or self.system.q_dep[:] == [q[3]]
        assert 'q' in exclude or self.system.q[:] == q[:4]
        assert 'u_ind' in exclude or self.system.u_ind[:] == u[:2]
        assert 'u_dep' in exclude or self.system.u_dep[:] == u[2:4]
        assert 'u' in exclude or self.system.u[:] == u[:4]
        assert 'kdes' in exclude or self.system.kdes[:] == [
            ui - qdi for ui, qdi in zip(u[:4], qd[:4])]
        assert ('holonomic_constraints' in exclude or
                self.system.holonomic_constraints[:] == [q[2] - q[0] + q[1]])
        assert ('nonholonomic_constraints' in exclude or
                self.system.nonholonomic_constraints[:] == [u[3] - qd[1] + u[2]]
                )
        assert ('bodies' in exclude or
                self.system.bodies == tuple(self.bodies))
        assert ('joints' in exclude or
                self.system.joints == tuple(self.joints))

    def test_filled_system(self, _filled_system_setup):
        self._filled_system_check()

    @pytest.mark.parametrize('origin', [None, Point('origin')])
    @pytest.mark.parametrize('frame', [None, ReferenceFrame('frame')])
    def test_init(self, origin, frame):
        if origin is None and frame is None:
            self.system = System()
        else:
            self.system = System(origin, frame)
        if origin is None:
            assert self.system.origin.name == 'inertial_origin'
        else:
            assert self.system.origin == origin
        if frame is None:
            assert self.system.frame.name == 'inertial_frame'
        else:
            assert self.system.frame == frame
        self._empty_system_check()
        assert isinstance(self.system.q_ind, ImmutableMatrix)
        assert isinstance(self.system.q_dep, ImmutableMatrix)
        assert isinstance(self.system.q, ImmutableMatrix)
        assert isinstance(self.system.u_ind, ImmutableMatrix)
        assert isinstance(self.system.u_dep, ImmutableMatrix)
        assert isinstance(self.system.u, ImmutableMatrix)
        assert isinstance(self.system.kdes, ImmutableMatrix)
        assert isinstance(self.system.holonomic_constraints, ImmutableMatrix)
        assert isinstance(self.system.nonholonomic_constraints, ImmutableMatrix)

    def test_from_newtonian_rigid_body(self):
        rb = RigidBody('body')
        self.system = System.from_newtonian(rb)
        assert self.system.origin == rb.masscenter
        assert self.system.frame == rb.frame
        self._empty_system_check(exclude=('bodies',))
        self.system.bodies = (rb,)

    def test_from_newtonian_particle(self):
        pt = Particle('particle')
        with pytest.raises(TypeError):
            System.from_newtonian(pt)

    @pytest.mark.parametrize('args, kwargs, exp_q_ind, exp_q_dep, exp_q', [
        (q[:3], {}, q[:3], [], q[:3]),
        (q[:3], {'independent': True}, q[:3], [], q[:3]),
        (q[:3], {'independent': False}, [], q[:3], q[:3]),
        (q[:3], {'independent': [True, False, True]}, [q[0], q[2]], [q[1]],
         [q[0], q[2], q[1]]),
    ])
    def test_coordinates(self, _empty_system_setup, args, kwargs,
                         exp_q_ind, exp_q_dep, exp_q):
        # Test add_coordinates
        self.system.add_coordinates(*args, **kwargs)
        assert self.system.q_ind[:] == exp_q_ind
        assert self.system.q_dep[:] == exp_q_dep
        assert self.system.q[:] == exp_q
        self._empty_system_check(exclude=('q_ind', 'q_dep', 'q'))
        # Test setter for q_ind and q_dep
        self.system.q_ind = exp_q_ind
        self.system.q_dep = exp_q_dep
        assert self.system.q_ind[:] == exp_q_ind
        assert self.system.q_dep[:] == exp_q_dep
        assert self.system.q[:] == exp_q
        self._empty_system_check(exclude=('q_ind', 'q_dep', 'q'))

    @pytest.mark.parametrize('func', ['add_coordinates', 'add_speeds'])
    @pytest.mark.parametrize('args, kwargs', [
        ((q[0], q[5]), {}),
        ((u[0], u[5]), {}),
        ((q[0],), {'independent': False}),
        ((u[0],), {'independent': False}),
        ((u[0], q[5]), {}),
        ((symbols('a'), q[5]), {}),
    ])
    def test_coordinates_speeds_invalid(self, _filled_system_setup, func, args,
                                        kwargs):
        with pytest.raises(ValueError):
            getattr(self.system, func)(*args, **kwargs)
        self._filled_system_check()

    @pytest.mark.parametrize('args, kwargs, exp_u_ind, exp_u_dep, exp_u', [
        (u[:3], {}, u[:3], [], u[:3]),
        (u[:3], {'independent': True}, u[:3], [], u[:3]),
        (u[:3], {'independent': False}, [], u[:3], u[:3]),
        (u[:3], {'independent': [True, False, True]}, [u[0], u[2]], [u[1]],
         [u[0], u[2], u[1]]),
    ])
    def test_speeds(self, _empty_system_setup, args, kwargs, exp_u_ind,
                    exp_u_dep, exp_u):
        # Test add_speeds
        self.system.add_speeds(*args, **kwargs)
        assert self.system.u_ind[:] == exp_u_ind
        assert self.system.u_dep[:] == exp_u_dep
        assert self.system.u[:] == exp_u
        self._empty_system_check(exclude=('u_ind', 'u_dep', 'u'))
        # Test setter for u_ind and u_dep
        self.system.u_ind = exp_u_ind
        self.system.u_dep = exp_u_dep
        assert self.system.u_ind[:] == exp_u_ind
        assert self.system.u_dep[:] == exp_u_dep
        assert self.system.u[:] == exp_u
        self._empty_system_check(exclude=('u_ind', 'u_dep', 'u'))

    @pytest.mark.parametrize('prop, add_func, args, kwargs', [
        ('q_ind', 'add_coordinates', (q[0],), {}),
        ('q_dep', 'add_coordinates', (q[3],), {'independent': False}),
        ('u_ind', 'add_speeds', (u[0],), {}),
        ('u_dep', 'add_speeds', (u[3],), {'independent': False}),
        ('kdes', 'add_kdes', (qd[0] - u[0],), {}),
        ('holonomic_constraints', 'add_holonomic_constraints',
         (q[0] - q[1],), {}),
        ('nonholonomic_constraints', 'add_nonholonomic_constraints',
         (u[0] - u[1],), {}),
        ('bodies', 'add_bodies', (RigidBody('body'),), {}),
        ('loads', 'add_loads', (Force(Point('P'), ReferenceFrame('N').x),), {}),
    ])
    def test_add_after_reset(self, _filled_system_setup, prop, add_func, args,
                             kwargs):
        setattr(self.system, prop, ())
        self._filled_system_check(exclude=(prop, 'q', 'u'))
        assert list(getattr(self.system, prop)[:]) == []
        getattr(self.system, add_func)(*args, **kwargs)
        assert list(getattr(self.system, prop)[:]) == list(args)

    @pytest.mark.parametrize('prop, add_func, value, error', [
        ('q_ind', 'add_coordinates', symbols('a'), ValueError),
        ('q_dep', 'add_coordinates', symbols('a'), ValueError),
        ('u_ind', 'add_speeds', symbols('a'), ValueError),
        ('u_dep', 'add_speeds', symbols('a'), ValueError),
        ('kdes', 'add_kdes', 7, TypeError),
        ('holonomic_constraints', 'add_holonomic_constraints', 7, TypeError),
        ('nonholonomic_constraints', 'add_nonholonomic_constraints', 7,
         TypeError),
        ('bodies', 'add_bodies', symbols('a'), TypeError),
        ('loads', 'add_loads', symbols('a'), TypeError),
    ])
    def test_type_error(self, _filled_system_setup, prop, add_func, value,
                        error):
        with pytest.raises(error):
            getattr(self.system, add_func)(value)
        with pytest.raises(error):
            setattr(self.system, prop, value)
        self._filled_system_check()

    @pytest.mark.parametrize('args, kwargs, exp_kdes', [
        ((), {}, [ui - qdi for ui, qdi in zip(u[:4], qd[:4])]),
        ((u[4] - qd[4], u[5] - qd[5]), {},
         [ui - qdi for ui, qdi in zip(u[:6], qd[:6])]),
    ])
    def test_kdes(self, _filled_system_setup, args, kwargs, exp_kdes):
        # Test add_speeds
        self.system.add_kdes(*args, **kwargs)
        self._filled_system_check(exclude=('kdes',))
        assert self.system.kdes[:] == exp_kdes
        # Test setter for kdes
        self.system.kdes = exp_kdes
        self._filled_system_check(exclude=('kdes',))
        assert self.system.kdes[:] == exp_kdes

    @pytest.mark.parametrize('args, kwargs', [
        ((u[0] - qd[0], u[4] - qd[4]), {}),
        ((-(u[0] - qd[0]), u[4] - qd[4]), {}),
        (([u[0] - u[0], u[4] - qd[4]]), {}),
    ])
    def test_kdes_invalid(self, _filled_system_setup, args, kwargs):
        with pytest.raises(ValueError):
            self.system.add_kdes(*args, **kwargs)
        self._filled_system_check()

    @pytest.mark.parametrize('args, kwargs, exp_con', [
        ((), {}, [q[2] - q[0] + q[1]]),
        ((q[4] - q[5], q[5] + q[3]), {},
         [q[2] - q[0] + q[1], q[4] - q[5], q[5] + q[3]]),
    ])
    def test_holonomic_constraints(self, _filled_system_setup, args, kwargs,
                                   exp_con):
        # Test add_holonomic_constraints
        self.system.add_holonomic_constraints(*args, **kwargs)
        self._filled_system_check(exclude=('holonomic_constraints',))
        assert self.system.holonomic_constraints[:] == exp_con
        # Test setter for holonomic_constraints
        self.system.holonomic_constraints = exp_con
        self._filled_system_check(exclude=('holonomic_constraints',))
        assert self.system.holonomic_constraints[:] == exp_con

    @pytest.mark.parametrize('args, kwargs', [
        ((q[2] - q[0] + q[1], q[4] - q[3]), {}),
        ((-(q[2] - q[0] + q[1]), q[4] - q[3]), {}),
        ((q[0] - q[0], q[4] - q[3]), {}),
    ])
    def test_holonomic_constraints_invalid(self, _filled_system_setup, args,
                                           kwargs):
        with pytest.raises(ValueError):
            self.system.add_holonomic_constraints(*args, **kwargs)
        self._filled_system_check()

    @pytest.mark.parametrize('args, kwargs, exp_con', [
        ((), {}, [u[3] - qd[1] + u[2]]),
        ((u[4] - u[5], u[5] + u[3]), {},
         [u[3] - qd[1] + u[2], u[4] - u[5], u[5] + u[3]]),
    ])
    def test_nonholonomic_constraints(self, _filled_system_setup, args, kwargs,
                                      exp_con):
        # Test add_nonholonomic_constraints
        self.system.add_nonholonomic_constraints(*args, **kwargs)
        self._filled_system_check(exclude=('nonholonomic_constraints',))
        assert self.system.nonholonomic_constraints[:] == exp_con
        # Test setter for nonholonomic_constraints
        self.system.nonholonomic_constraints = exp_con
        self._filled_system_check(exclude=('nonholonomic_constraints',))
        assert self.system.nonholonomic_constraints[:] == exp_con

    @pytest.mark.parametrize('args, kwargs', [
        ((u[3] - qd[1] + u[2], u[4] - u[3]), {}),
        ((-(u[3] - qd[1] + u[2]), u[4] - u[3]), {}),
        ((u[0] - u[0], u[4] - u[3]), {}),
        (([u[0] - u[0], u[4] - u[3]]), {}),
    ])
    def test_nonholonomic_constraints_invalid(self, _filled_system_setup, args,
                                              kwargs):
        with pytest.raises(ValueError):
            self.system.add_nonholonomic_constraints(*args, **kwargs)
        self._filled_system_check()

    def test_bodies(self, _filled_system_setup):
        rb1, rb2 = RigidBody('rb1'), RigidBody('rb2')
        p1, p2 = Particle('p1'), Particle('p2')
        self.system.add_bodies(rb1, p1)
        assert self.system.bodies == (*self.bodies, rb1, p1)
        self.system.add_bodies(p2)
        assert self.system.bodies == (*self.bodies, rb1, p1, p2)
        self.system.bodies = []
        assert self.system.bodies == ()
        self.system.bodies = p2
        assert self.system.bodies == (p2,)
        symb = symbols('symb')
        pytest.raises(TypeError, lambda: self.system.add_bodies(symb))
        pytest.raises(ValueError, lambda: self.system.add_bodies(p2))
        with pytest.raises(TypeError):
            self.system.bodies = (rb1, rb2, p1, p2, symb)
        assert self.system.bodies == (p2,)

    def test_add_loads(self):
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
        assert system.remove_load(rb1.masscenter) == ((rb1.masscenter, N.x),)
        assert system.remove_load(rb2.masscenter) == ((rb2.masscenter, -N.x),)
        assert system.remove_load(rb2.masscenter) == ()
        system.apply_torque(rb1, N.x, rb2)
        assert system.loads == ((A, A.x), (N, N.x), (A, -N.x))
        system.apply_force(p2, A.x)
        assert system.loads == ((A, A.x), (N, N.x), (A, -N.x), (mc2, A.x))
        pytest.raises(ValueError, lambda: system.add_loads((N, N.x, N.y)))
        with pytest.raises(TypeError):
            system.loads = (N, N.x)
        assert system.loads == ((A, A.x), (N, N.x), (A, -N.x), (mc2, A.x))
        system.clear_loads()
        assert system.loads == ()

    def test_remove_load(self):
        system = System()
        N = ReferenceFrame('N')
        pnt = Point('pnt')
        rb = RigidBody('rb')
        p = Particle('p')
        system.add_loads((rb.masscenter, N.x), (rb.frame, N.x), (pnt, N.z),
                         (p.masscenter, N.z), (pnt, N.x), (N, N.x), (N, N.y),
                         (N, N.z))
        assert system.loads == (
        (rb.masscenter, N.x), (rb.frame, N.x), (pnt, N.z),
        (p.masscenter, N.z), (pnt, N.x), (N, N.x), (N, N.y),
        (N, N.z))
        assert system.remove_load(rb) == ((rb.masscenter, N.x), (rb.frame, N.x))
        assert system.loads == ((pnt, N.z), (p.masscenter, N.z), (pnt, N.x),
                                (N, N.x), (N, N.y), (N, N.z))
        assert system.remove_load(p) == ((p.masscenter, N.z),)
        assert system.loads == ((pnt, N.z), (pnt, N.x), (N, N.x), (N, N.y),
                                (N, N.z))
        assert system.remove_load(pnt) == ((pnt, N.z), (pnt, N.x))
        assert system.loads == ((N, N.x), (N, N.y), (N, N.z))
        assert system.remove_load(pnt) == ()
        assert system.remove_load(N) == ((N, N.x), (N, N.y), (N, N.z))
        assert system.loads == ()
        assert system.remove_load(N) == ()

    def test_add_joints(self):
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
        assert system.q_ind == ImmutableMatrix([q1])
        assert system.u_ind == ImmutableMatrix([u1])
        assert system.kdes == ImmutableMatrix([u1 - q1.diff(t)])
        system.add_bodies(rb4)
        system.add_coordinates(q3)
        system.add_kdes(u3 - q3.diff(t))
        system.add_joints(J3)
        assert system.joints == (J1, J3)
        assert system.bodies == (rb1, rb2, rb4, rb3)
        assert system.q_ind == ImmutableMatrix([q1, q3])
        assert system.u_ind == ImmutableMatrix([u1, u3])
        assert system.kdes == ImmutableMatrix(
            [u1 - q1.diff(t), u3 - q3.diff(t)])
        system.add_kdes(-(u2 - q2.diff(t)))
        system.add_joints(J2)
        assert system.joints == (J1, J3, J2)
        assert system.bodies == (rb1, rb2, rb4, rb3)
        assert system.q_ind == ImmutableMatrix([q1, q3, q2])
        assert system.u_ind == ImmutableMatrix([u1, u3, u2])
        assert system.kdes == ImmutableMatrix([u1 - q1.diff(t), u3 - q3.diff(t),
                                               -(u2 - q2.diff(t))])
        system.add_joints(J_lag)
        assert system.joints == (J1, J3, J2, J_lag)
        assert system.bodies == (rb1, rb2, rb4, rb3, rb5)
        assert system.q_ind == ImmutableMatrix([q1, q3, q2, q4])
        assert system.u_ind == ImmutableMatrix([u1, u3, u2, q4.diff(t)])
        assert system.kdes == ImmutableMatrix([u1 - q1.diff(t), u3 - q3.diff(t),
                                               -(u2 - q2.diff(t))])
        assert system.q_dep[:] == []
        assert system.u_dep[:] == []
        pytest.raises(ValueError, lambda: system.add_joints(J2))
        pytest.raises(TypeError, lambda: system.add_joints(rb1))

    @pytest.mark.parametrize('name, joint_index', [
        ('J1', 0),
        ('J2', 1),
        ('not_existing', None),
    ])
    def test_get_joint(self, _filled_system_setup, name, joint_index):
        joint = self.system.get_joint(name)
        if joint_index is None:
            assert joint is None
        else:
            assert joint == self.joints[joint_index]

    @pytest.mark.parametrize('name, body_index', [
        ('rb1', 0),
        ('rb3', 2),
        ('not_existing', None),
    ])
    def test_get_body(self, _filled_system_setup, name, body_index):
        body = self.system.get_body(name)
        if body_index is None:
            assert body is None
        else:
            assert body == self.bodies[body_index]
