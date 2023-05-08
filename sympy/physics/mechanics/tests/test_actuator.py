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
    LinearSpring,
    TorqueActuator,
)
from sympy.physics.mechanics._pathway import LinearPathway, PathwayBase

if USE_SYMENGINE:
    from sympy.core.backend import Basic as ExprType
else:
    from sympy.core.expr import Expr as ExprType


target = RigidBody('target')
reaction = RigidBody('reaction')


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


class TestLinearSpring:

    @pytest.fixture(autouse=True)
    def _linear_pathway_fixture(self) -> None:
        self.stiffness = Symbol('k')
        self.l = Symbol('l')
        self.pA = Point('pA')
        self.pB = Point('pB')
        self.pathway = LinearPathway(self.pA, self.pB)
        self.q = dynamicsymbols('q')
        self.N = ReferenceFrame('N')

    def test_is_force_actuator_subclass(self) -> None:
        assert issubclass(LinearSpring, ForceActuator)

    def test_is_actuator_base_subclass(self) -> None:
        assert issubclass(LinearSpring, ActuatorBase)

    @pytest.mark.parametrize(
        'stiffness, equilibrium_length, force',
        [
            (Symbol('k'), S.Zero, -Symbol('k')*sqrt(dynamicsymbols('q')**2)),
            (Symbol('k'), Symbol('l'), -Symbol('k')*(sqrt(dynamicsymbols('q')**2) - Symbol('l'))),
        ]
    )
    def test_valid_constructor(
        self,
        stiffness: ExprType,
        equilibrium_length: ExprType,
        force: ExprType,
    ) -> None:
        self.pB.set_pos(self.pA, self.q * self.N.x)
        spring = LinearSpring(
            stiffness,
            self.pathway,
            equilibrium_length=equilibrium_length,
        )

        assert isinstance(spring, LinearSpring)

        assert hasattr(spring, 'stiffness')
        assert isinstance(spring.stiffness, ExprType)
        assert spring.stiffness == stiffness

        assert hasattr(spring, 'pathway')
        assert isinstance(spring.pathway, LinearPathway)
        assert spring.pathway == self.pathway

        assert hasattr(spring, 'equilibrium_length')
        assert isinstance(spring.equilibrium_length, ExprType)
        assert spring.equilibrium_length == equilibrium_length

        assert hasattr(spring, 'force')
        assert isinstance(spring.force, ExprType)
        assert spring.force == force


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
        self.target = RigidBody('target', frame=self.N)
        self.reaction = RigidBody('reaction', frame=self.A)

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
        'target_frame, reaction_frame',
        [
            (target.frame, reaction.frame),
            (target, reaction.frame),
            (target.frame, reaction),
            (target, reaction),
        ]
    )
    def test_valid_constructor_with_reaction(
        self,
        torque: ExprType,
        target_frame: ReferenceFrame | RigidBody,
        reaction_frame: ReferenceFrame | RigidBody,
    ) -> None:
        instance = TorqueActuator(torque, self.N.z, target_frame, reaction_frame)
        assert isinstance(instance, TorqueActuator)

        assert hasattr(instance, 'torque')
        assert isinstance(instance.torque, ExprType)
        assert instance.torque == torque

        assert hasattr(instance, 'axis')
        assert isinstance(instance.axis, Vector)
        assert instance.axis == self.N.z

        assert hasattr(instance, 'target_frame')
        assert isinstance(instance.target_frame, ReferenceFrame)
        assert instance.target_frame == target.frame

        assert hasattr(instance, 'reaction_frame')
        assert isinstance(instance.reaction_frame, ReferenceFrame)
        assert instance.reaction_frame == reaction.frame

    @pytest.mark.parametrize(
        'torque',
        [
            Symbol('T'),
            dynamicsymbols('T'),
            Symbol('T')**2 + Symbol('T'),
        ]
    )
    @pytest.mark.parametrize('target_frame', [target.frame, target])
    def test_valid_constructor_without_reaction(
        self,
        torque: ExprType,
        target_frame: ReferenceFrame | RigidBody,
    ) -> None:
        instance = TorqueActuator(torque, self.N.z, target_frame)
        assert isinstance(instance, TorqueActuator)

        assert hasattr(instance, 'torque')
        assert isinstance(instance.torque, ExprType)
        assert instance.torque == torque

        assert hasattr(instance, 'axis')
        assert isinstance(instance.axis, Vector)
        assert instance.axis == self.N.z

        assert hasattr(instance, 'target_frame')
        assert isinstance(instance.target_frame, ReferenceFrame)
        assert instance.target_frame == target.frame

        assert hasattr(instance, 'reaction_frame')
        assert instance.reaction_frame is None

    @pytest.mark.parametrize(
        'axis',
        [
            Symbol('a'),
            dynamicsymbols('a'),
        ]
    )
    def test_invalid_constructor_axis_not_vector(self, axis: Any) -> None:
        with pytest.raises(TypeError):
            _ = TorqueActuator(self.torque, axis, self.target, self.reaction)  # type: ignore

    @pytest.mark.parametrize(
        'frames',
        [
            (None, ReferenceFrame('child')),
            (ReferenceFrame('parent'), True),
            (None, RigidBody('child')),
            (RigidBody('parent'), True),
        ]
    )
    def test_invalid_frames_not_frame(self, frames: Sequence[Any]) -> None:
        with pytest.raises(TypeError):
            _ = TorqueActuator(self.torque, self.N.z, *frames)  # type: ignore

    def test_repr_without_reaction(self) -> None:
        actuator = TorqueActuator(self.torque, self.N.z, self.target)
        expected = 'TorqueActuator(T, axis=N.z, target_frame=N)'
        assert repr(actuator) == expected

    def test_repr_with_reaction(self) -> None:
        actuator = TorqueActuator(self.torque, self.N.z, self.target, self.reaction)
        expected = 'TorqueActuator(T, axis=N.z, target_frame=N, reaction_frame=A)'
        assert repr(actuator) == expected

    def test_at_pin_joint_constructor(self) -> None:
        pin_joint = PinJoint(
            'pin',
            self.target,
            self.reaction,
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

        assert hasattr(instance, 'target_frame')
        assert isinstance(instance.target_frame, ReferenceFrame)
        assert instance.target_frame == self.N

        assert hasattr(instance, 'reaction_frame')
        assert isinstance(instance.reaction_frame, ReferenceFrame)
        assert instance.reaction_frame == self.A

    def test_at_pin_joint_pin_joint_not_pin_joint_invalid(self) -> None:
        with pytest.raises(TypeError):
            _ = TorqueActuator.at_pin_joint(self.torque, Symbol('pin'))  # type: ignore

    def test_to_loads_without_reaction(self) -> None:
        actuator = TorqueActuator(self.torque, self.N.z, self.target)
        expected = [
            (self.N, self.torque * self.N.z),
        ]
        assert actuator.to_loads() == expected

    def test_to_loads_with_reaction(self) -> None:
        actuator = TorqueActuator(self.torque, self.N.z, self.target, self.reaction)
        expected = [
            (self.N, self.torque * self.N.z),
            (self.A, - self.torque * self.N.z),
        ]
        assert actuator.to_loads() == expected
