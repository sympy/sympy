import pytest

from sympy.core.backend import ImmutableMatrix, symbols
from sympy.physics.mechanics import (
    System, ReferenceFrame, Point, RigidBody, Particle, dynamicsymbols,
    PinJoint, PrismaticJoint)

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
        assert 'q_ind' not in exclude or self.system.q_ind[:] == q[:3]
        assert 'q_dep' not in exclude or self.system.q_dep[:] == [q[3]]
        assert 'q' not in exclude or self.system.q[:] == q[:4]
        assert 'u_ind' not in exclude or self.system.u_ind[:] == u[:2]
        assert 'u_dep' not in exclude or self.system.u_dep[:] == u[2:4]
        assert 'u' not in exclude or self.system.u[:] == u[:4]
        assert 'kdes' not in exclude or self.system.kdes[:] == [
            ui - qdi for ui, qdi in zip(u[:4], qd[:4])]
        assert ('holonomic_constraints' not in exclude or
                self.system.holonomic_constraints[:] == [q[2] - q[0] + q[1]])
        assert ('nonholonomic_constraints' not in exclude or
                self.system.nonholonomic_constraints[:] == [u[3] - qd[1] + u[2]]
                )
        assert ('bodies' not in exclude or
                self.system.bodies == tuple(self.bodies))
        assert ('joints' not in exclude or
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
    def test_coordinates_simple(self, _empty_system_setup, args, kwargs,
                                exp_q_ind, exp_q_dep, exp_q):
        self._empty_system_check()
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

    def test_coordinate_additional(self, _filled_system_setup):
        self.system.add_coordinates(q[4], independent=False)
        assert self.system.q_ind[:] == q[:3]
        assert self.system.q_dep[:] == q[3:5]
        assert self.system.q[:] == q[:5]
        # Test resetting the coordinates
        self.system.q_ind = []
        assert self.system.q_ind[:] == []
        assert self.system.q_dep[:] == q[3:5]
        assert self.system.q[:] == q[3:5]
        self.system.q_dep = []
        assert self.system.q_ind[:] == []
        assert self.system.q_dep[:] == []
        assert self.system.q[:] == []
        # Test if filling the coordinates after reset works
        self.system.add_coordinates(q[0])
        assert self.system.q_ind[:] == [q[0]]
        assert self.system.q_dep[:] == []
        assert self.system.q[:] == [q[0]]
        self.system.add_coordinates(q[1], independent=False)
        assert self.system.q_ind[:] == [q[0]]
        assert self.system.q_dep[:] == [q[1]]
        assert self.system.q[:] == q[:2]

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
    def test_speeds_simple(self, _empty_system_setup, args, kwargs, exp_u_ind,
                           exp_u_dep, exp_u):
        self._empty_system_check()
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

    def test_speed_additional(self, _filled_system_setup):
        self.system.add_speeds(u[4], independent=False)
        assert self.system.u_ind[:] == u[:2]
        assert self.system.u_dep[:] == u[2:5]
        assert self.system.u[:] == u[:5]
        # Test resetting the coordinates
        self.system.u_ind = []
        assert self.system.u_ind[:] == []
        assert self.system.u_dep[:] == u[2:5]
        assert self.system.u[:] == u[2:5]
        self.system.u_dep = []
        assert self.system.u_ind[:] == []
        assert self.system.u_dep[:] == []
        assert self.system.u[:] == []
        # Test if filling the coordinates after reset works
        self.system.add_speeds(u[0])
        assert self.system.u_ind[:] == [u[0]]
        assert self.system.u_dep[:] == []
        assert self.system.u[:] == [u[0]]
        self.system.add_speeds(u[1], independent=False)
        assert self.system.u_ind[:] == [u[0]]
        assert self.system.u_dep[:] == [u[1]]
        assert self.system.u[:] == u[:2]
