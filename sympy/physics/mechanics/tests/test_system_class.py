import pytest

from sympy.core.backend import ImmutableMatrix
from sympy.physics.mechanics import (
    System, ReferenceFrame, Point, RigidBody, Particle, dynamicsymbols)

t = dynamicsymbols._t  # type: ignore
q = dynamicsymbols('q:6')  # type: ignore
qd = dynamicsymbols('q:6', 1)  # type: ignore
u = dynamicsymbols('u:6')  # type: ignore


class TestSystem:
    @staticmethod
    def _empty_system(system, exclude=()):
        matrices = ('q_ind', 'q_dep', 'q', 'u_ind', 'u_dep', 'u', 'kdes',
                    'holonomic_constraints', 'nonholonomic_constraints')
        tuples = ('loads', 'bodies', 'joints')
        for attr in matrices:
            if attr not in exclude:
                assert getattr(system, attr)[:] == []
        for attr in tuples:
            if attr not in exclude:
                assert getattr(system, attr) == ()
        if 'eom_method' not in exclude:
            assert system.eom_method is None

    @pytest.mark.parametrize('origin', [None, Point('origin')])
    @pytest.mark.parametrize('frame', [None, ReferenceFrame('frame')])
    def test_init(self, origin, frame):
        if origin is None and frame is None:
            system = System()
        else:
            system = System(origin, frame)
        if origin is None:
            assert system.origin.name == 'inertial_origin'
        else:
            assert system.origin == origin
        if frame is None:
            assert system.frame.name == 'inertial_frame'
        else:
            assert system.frame == frame
        self._empty_system(system)
        assert isinstance(system.q_ind, ImmutableMatrix)
        assert isinstance(system.q_dep, ImmutableMatrix)
        assert isinstance(system.q, ImmutableMatrix)
        assert isinstance(system.u_ind, ImmutableMatrix)
        assert isinstance(system.u_dep, ImmutableMatrix)
        assert isinstance(system.u, ImmutableMatrix)
        assert isinstance(system.kdes, ImmutableMatrix)
        assert isinstance(system.holonomic_constraints, ImmutableMatrix)
        assert isinstance(system.nonholonomic_constraints, ImmutableMatrix)

    def test_from_newtonian_rigid_body(self):
        rb = RigidBody('body')
        system = System.from_newtonian(rb)
        assert system.origin == rb.masscenter
        assert system.frame == rb.frame
        self._empty_system(system, exclude=('bodies',))
        system.bodies = (rb,)

    def test_from_newtonian_particle(self):
        pt = Particle('particle')
        with pytest.raises(TypeError):
            System.from_newtonian(pt)

    @pytest.fixture(autouse=True)
    def _empty_setup(self):
        self.system = System(Point('origin'), ReferenceFrame('frame'))

    @pytest.mark.parametrize('args, kwargs, exp_q_ind, exp_q_dep, exp_q', [
        (q[:3], {}, q[:3], [], q[:3]),
        (q[:3], {'independent': True}, q[:3], [], q[:3]),
        (q[:3], {'independent': False}, [], q[:3], q[:3]),
        (q[:3], {'independent': [True, False, True]}, [q[0], q[2]], [q[1]],
         [q[0], q[2], q[1]]),
    ])
    def test_coordinates_simple(self, _empty_setup, args, kwargs, exp_q_ind,
                                exp_q_dep, exp_q):
        self._empty_system(self.system)
        # Test add_coordinates
        self.system.add_coordinates(*args, **kwargs)
        assert self.system.q_ind[:] == exp_q_ind
        assert self.system.q_dep[:] == exp_q_dep
        assert self.system.q[:] == exp_q
        self._empty_system(self.system, exclude=('q_ind', 'q_dep', 'q'))
        # Test setter for q_ind and q_dep
        self.system.q_ind = exp_q_ind
        self.system.q_dep = exp_q_dep
        assert self.system.q_ind[:] == exp_q_ind
        assert self.system.q_dep[:] == exp_q_dep
        assert self.system.q[:] == exp_q
        self._empty_system(self.system, exclude=('q_ind', 'q_dep', 'q'))

    @pytest.mark.parametrize('args, kwargs, exp_u_ind, exp_u_dep, exp_u', [
        (u[:3], {}, u[:3], [], u[:3]),
        (u[:3], {'independent': True}, u[:3], [], u[:3]),
        (u[:3], {'independent': False}, [], u[:3], u[:3]),
        (u[:3], {'independent': [True, False, True]}, [u[0], u[2]], [u[1]],
         [u[0], u[2], u[1]]),
    ])
    def test_speeds_simple(self, _empty_setup, args, kwargs, exp_u_ind,
                           exp_u_dep, exp_u):
        self._empty_system(self.system)
        # Test add_speeds
        self.system.add_speeds(*args, **kwargs)
        assert self.system.u_ind[:] == exp_u_ind
        assert self.system.u_dep[:] == exp_u_dep
        assert self.system.u[:] == exp_u
        self._empty_system(self.system, exclude=('u_ind', 'u_dep', 'u'))
        # Test setter for u_ind and u_dep
        self.system.u_ind = exp_u_ind
        self.system.u_dep = exp_u_dep
        assert self.system.u_ind[:] == exp_u_ind
        assert self.system.u_dep[:] == exp_u_dep
        assert self.system.u[:] == exp_u
        self._empty_system(self.system, exclude=('u_ind', 'u_dep', 'u'))
