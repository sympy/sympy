import pytest
from sympy.physics.mechanics import (
    System, ReferenceFrame, Point, RigidBody, Particle)


class TestSystem:
    @staticmethod
    def _empty_system(system, exclude=()):
        matrices = ("q_ind", "q_dep", "q", "u_ind", "u_dep", "u", "kdes",
                    "holonomic_constraints", "nonholonomic_constraints")
        tuples = ("loads", "bodies", "joints")
        for attr in matrices:
            if attr not in exclude:
                assert getattr(system, attr)[:] == []
        for attr in tuples:
            if attr not in exclude:
                assert getattr(system, attr) == ()
        if 'eom_method' not in exclude:
            assert system.eom_method is None

    @pytest.mark.parametrize("origin", [None, Point("origin")])
    @pytest.mark.parametrize("frame", [None, ReferenceFrame('frame')])
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

