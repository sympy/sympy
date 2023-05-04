"""Tests for the ``sympy.physics.mechanics._actuator.py`` module."""

from __future__ import annotations

from typing import Any, Sequence

import pytest

from sympy.core.backend import USE_SYMENGINE, Matrix, Symbol, sqrt
from sympy.physics.mechanics import (
    KanesMethod,
    Particle,
    PinJoint,
    Point,
    ReferenceFrame,
    RigidBody,
    Vector,
    dynamicsymbols,
)
from sympy.physics.mechanics._actuator import (
    ActuatorBase,
    ForceActuator,
    TorqueActuator,
)
from sympy.physics.mechanics._pathway import LinearPathway, PathwayBase

if USE_SYMENGINE:
    from sympy.core.backend import Basic as ExprType
else:
    from sympy.core.expr import Expr as ExprType


parent = RigidBody('parent')
child = RigidBody('child')


class TestForceActuator:

    @pytest.fixture(autouse=True)
    def _linear_pathway_fixture(self) -> None:
        self.force = Symbol('F')
        self.pA = Point('pA')
        self.pB = Point('pB')
        self.pathway = LinearPathway(self.pA, self.pB)
        self.q1 = dynamicsymbols('q1')
        self.q2 = dynamicsymbols('q2')
        self.q3 = dynamicsymbols('q3')
        self.q1d = dynamicsymbols('q1', 1)
        self.q2d = dynamicsymbols('q2', 1)
        self.q3d = dynamicsymbols('q3', 1)
        self.N = ReferenceFrame('N')

    def test_is_actuator_base_subclass(self) -> None:
        assert issubclass(ForceActuator, ActuatorBase)

    @pytest.mark.parametrize(
        'force',
        [
            Symbol('F'),
            dynamicsymbols('F'),
            Symbol('F')**2 + Symbol('F'),
        ]
    )
    def test_valid_constructor_force(self, force: ExprType) -> None:
        instance = ForceActuator(force, self.pathway)
        assert isinstance(instance, ForceActuator)
        assert hasattr(instance, 'force')
        assert isinstance(instance.force, ExprType)
        assert instance.force == force

    def test_invalid_constructor_force_not_expr(self) -> None:
        with pytest.raises(TypeError):
            _ = ForceActuator(None, self.pathway)  # type: ignore

    @pytest.mark.parametrize(
        'pathway',
        [
            LinearPathway(Point('pA'), Point('pB')),
        ]
    )
    def test_valid_constructor_pathway(self, pathway: PathwayBase) -> None:
        instance = ForceActuator(self.force, pathway)
        assert isinstance(instance, ForceActuator)
        assert hasattr(instance, 'pathway')
        assert isinstance(instance.pathway, LinearPathway)
        assert instance.pathway == pathway

    def test_invalid_constructor_pathway_not_pathway_base(self) -> None:
        with pytest.raises(TypeError):
            _ = ForceActuator(self.force, None)  # type: ignore

    def test_repr(self) -> None:
        actuator = ForceActuator(self.force, self.pathway)
        expected = "ForceActuator(F, LinearPathway(pA, pB))"
        assert repr(actuator) == expected

    def test_to_loads_static_pathway(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.N.x)
        actuator = ForceActuator(self.force, self.pathway)
        expected = [
            (self.pA, - self.force * self.N.x),
            (self.pB, self.force * self.N.x),
        ]
        assert actuator.to_loads() == expected

    def test_to_loads_2D_pathway(self) -> None:
        self.pB.set_pos(self.pA, 2 * self.q1 * self.N.x)
        actuator = ForceActuator(self.force, self.pathway)
        expected = [
            (self.pA, - self.force * (self.q1 / sqrt(self.q1**2)) * self.N.x),
            (self.pB, self.force * (self.q1 / sqrt(self.q1**2)) * self.N.x),
        ]
        assert actuator.to_loads() == expected

    def test_to_loads_3D_pathway(self) -> None:
        self.pB.set_pos(
            self.pA,
            self.q1*self.N.x - self.q2*self.N.y + 2*self.q3*self.N.z,
        )
        actuator = ForceActuator(self.force, self.pathway)
        length = sqrt(self.q1**2 + self.q2**2 + 4*self.q3**2)
        pO_force = (
            - self.force * self.q1 * self.N.x / length
            + self.force * self.q2 * self.N.y / length
            - 2 * self.force * self.q3 * self.N.z / length
        )
        pI_force = (
            self.force * self.q1 * self.N.x / length
            - self.force * self.q2 * self.N.y / length
            + 2 * self.force * self.q3 * self.N.z / length
        )
        expected = [
            (self.pA, pO_force),
            (self.pB, pI_force),
        ]
        assert actuator.to_loads() == expected


def test_forced_mass_spring_damper_model():
    r"""A single degree of freedom translational forced mass-spring-damper.

    Notes
    =====

    This system is well known to have the governing equation:

    .. math::
        m \ddot{x} = F - k x - c \dot{x}

    where $F$ is an externally applied force, $m$ is the mass of the particle
    to which the spring and damper are attached, $k$ is the spring's stiffness,
    $c$ is the dampers damping coefficient, and $x$ is the generalized
    coordinate representing the system's single (translational) degree of
    freedom.

    """
    m = Symbol('m')
    k = Symbol('k')
    c = Symbol('c')
    F = Symbol('F')

    x = dynamicsymbols('x')
    dx = dynamicsymbols('x', 1)
    x_dot = dynamicsymbols('x_dot')

    frame = ReferenceFrame('N')
    origin = Point('pO')
    origin.set_vel(frame, 0)

    attachment = Point('pA')
    attachment.set_pos(origin, x * frame.x)

    mass = Particle('mass', attachment, m)
    pathway = LinearPathway(origin, attachment)
    stiffness = -k * pathway.length  # positive as assumes expansile force
    spring = ForceActuator(stiffness, pathway)
    damping = -c * pathway.extension_velocity  # negative as acts against velocity
    damper = ForceActuator(damping, pathway)

    kanes_method = KanesMethod(
        frame,
        q_ind=[x],
        u_ind=[x_dot],
        kd_eqs=[dx - x_dot],
    )
    bodies = [mass]
    loads = [
        (attachment, F * frame.x),
        *spring.to_loads(),
        *damper.to_loads(),
    ]
    kanes_method.kanes_equations(bodies, loads)

    assert kanes_method.mass_matrix == Matrix([[m]])
    assert kanes_method.forcing == Matrix([[F - c*x_dot - k*x]])


class TestTorqueActuator:

    @pytest.fixture(autouse=True)
    def _torque_actuator_fixture(self) -> None:
        self.torque = Symbol('T')
        self.N = ReferenceFrame('N')
        self.A = ReferenceFrame('A')
        self.parent = RigidBody('parent', frame=self.N)
        self.child = RigidBody('child', frame=self.A)

    def test_is_actuator_base_subclass(self) -> None:
        assert issubclass(TorqueActuator, ActuatorBase)

    @pytest.mark.parametrize(
        'torque',
        [
            Symbol('T'),
            dynamicsymbols('T'),
            Symbol('T')**2 + Symbol('T'),
        ]
    )
    @pytest.mark.parametrize(
        'frames',
        [
            (parent.frame, child.frame),
            (parent, child.frame),
            (parent.frame, child),
            (parent, child),
        ]
    )
    def test_valid_constructor(
        self,
        torque: ExprType,
        frames: tuple[ReferenceFrame | RigidBody, ReferenceFrame | RigidBody],
    ) -> None:
        instance = TorqueActuator(torque, self.N.z, frames[0], frames[1])
        assert isinstance(instance, TorqueActuator)

        assert hasattr(instance, 'torque')
        assert isinstance(instance.torque, ExprType)
        assert instance.torque == torque

        assert hasattr(instance, 'axis')
        assert isinstance(instance.axis, Vector)
        assert instance.axis == self.N.z

        assert hasattr(instance, 'frames')
        assert isinstance(instance.frames, tuple)
        assert isinstance(instance.frames[0], ReferenceFrame)
        assert isinstance(instance.frames[1], ReferenceFrame)
        assert instance.frames == (parent.frame, child.frame)

    @pytest.mark.parametrize(
        'axis',
        [
            Symbol('a'),
            dynamicsymbols('a'),
        ]
    )
    def test_invalid_constructor_axis_not_vector(self, axis: Any) -> None:
        with pytest.raises(TypeError):
            _ = TorqueActuator(self.torque, axis, self.parent, self.child)  # type: ignore

    @pytest.mark.parametrize(
        'frames',
        [
            (ReferenceFrame('A1'), ),
            (ReferenceFrame('A1'), ReferenceFrame('A2'), ReferenceFrame('A3')),
            (RigidBody('B1'), ),
            (RigidBody('B1'), RigidBody('B2'), RigidBody('B3')),
        ]
    )
    def test_invalid_frames_incorrect_number(
        self,
        frames: tuple[ReferenceFrame | RigidBody, ...],
    ) -> None:
        with pytest.raises(ValueError):
            _ = TorqueActuator(self.torque, self.N.z, *frames)  # type: ignore

    @pytest.mark.parametrize(
        'frames',
        [
            (None, ReferenceFrame('child')),
            (ReferenceFrame('parent'), None),
            (None, RigidBody('child')),
            (RigidBody('parent'), None),
        ]
    )
    def test_invalid_frames_not_frame(self, frames: Sequence[Any]) -> None:
        with pytest.raises(TypeError):
            _ = TorqueActuator(self.torque, self.N.z, *frames)  # type: ignore

    def test_repr(self) -> None:
        actuator = TorqueActuator(self.torque, self.N.z, self.parent, self.child)
        expected = "TorqueActuator(T, N.z, N, A)"
        assert repr(actuator) == expected

    def test_at_pin_joint_constructor(self) -> None:
        pin_joint = PinJoint(
            'pin',
            self.parent,
            self.child,
            coordinates=dynamicsymbols('q'),
            speeds=dynamicsymbols('u'),
            parent_interframe=self.N,
            joint_axis=self.N.z,
        )
        instance = TorqueActuator.at_pin_joint(self.torque, pin_joint)
        assert isinstance(instance, TorqueActuator)

        assert hasattr(instance, 'torque')
        assert isinstance(instance.torque, ExprType)
        assert instance.torque == self.torque

        assert hasattr(instance, 'axis')
        assert isinstance(instance.axis, Vector)
        assert instance.axis == self.N.z

        assert hasattr(instance, 'frames')
        assert isinstance(instance.frames, tuple)
        assert isinstance(instance.frames[0], ReferenceFrame)
        assert isinstance(instance.frames[1], ReferenceFrame)
        assert instance.frames == (self.N, self.A)

    def test_at_pin_joint_pin_joint_not_pin_joint_invalid(self) -> None:
        with pytest.raises(TypeError):
            _ = TorqueActuator.at_pin_joint(self.torque, Symbol('pin'))  # type: ignore
