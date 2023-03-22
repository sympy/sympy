import pytest

from sympy.core.backend import ImmutableMatrix, symbols
from sympy.physics.mechanics import (
    System, ReferenceFrame, Point, RigidBody, Particle, dynamicsymbols,
    PinJoint, PrismaticJoint, Force)

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
        (([q[0] - q[0], q[4] - q[3]]), {}),
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
        (([u[0] - u[0], u[4] - u[3]]), {}),
    ])
    def test_nonholonomic_constraints_invalid(self, _filled_system_setup, args,
                                              kwargs):
        with pytest.raises(ValueError):
            self.system.add_nonholonomic_constraints(*args, **kwargs)
        self._filled_system_check()
