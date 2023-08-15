"""Implementations of actuators for linked force and torque application.

Notes
=====

This module is experimental and so is named with a leading underscore to
indicate that the API is not yet stabilized and could be subject to breaking
changes.

"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

from sympy.core.backend import S, USE_SYMENGINE, sympify
from sympy.physics.mechanics import (
    PathwayBase,
    PinJoint,
    ReferenceFrame,
    RigidBody,
    Torque,
    Vector,
)

if USE_SYMENGINE:
    from sympy.core.backend import Basic as ExprType
else:
    from sympy.core.expr import Expr as ExprType

if TYPE_CHECKING:
    from sympy.physics.mechanics.loads import LoadBase


__all__ = [
    'ForceActuator',
    'LinearDamper',
    'LinearSpring',
    'TorqueActuator',
]


class ActuatorBase(ABC):
    """Abstract base class for all actuator classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom actuator types through subclassing.

    """

    def __init__(self) -> None:
        """Initializer for ``ActuatorBase``."""
        pass

    @abstractmethod
    def to_loads(self) -> list[LoadBase]:
        """Loads required by the equations of motion method classes.

        Explanation
        ===========

        ``KanesMethod`` requires a list of ``Point``-``Vector`` tuples to be
        passed to the ``loads`` parameters of its ``kanes_equations`` method
        when constructing the equations of motion. This method acts as a
        utility to produce the correctly-structred pairs of points and vectors
        required so that these can be easily concatenated with other items in
        the list of loads and passed to ``KanesMethod.kanes_equations``. These
        loads are also in the correct form to also be passed to the other
        equations of motion method classes, e.g. ``LagrangesMethod``.

        """
        pass

    def __repr__(self) -> str:
        """Default representation of an actuator."""
        return f'{self.__class__.__name__}()'


class ForceActuator(ActuatorBase):
    """Force-producing actuator.

    Explanation
    ===========

    A ``ForceActuator`` is an actuator that produces a (expansile) force along
    its length.

    Examples
    ========

    As the ``_actuator.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``ForceActuator`` must therefore be
    imported directly from the ``sympy.physics.mechanics._actuator`` module.

    >>> from sympy.physics.mechanics._actuator import ForceActuator

    To construct an actuator, an expression (or symbol) must be supplied to
    represent the force it can produce, alongside a pathway specifying its line
    of action. Let's also create a global reference frame and spatially fix one
    of the points in it while setting the other to be positioned such that it
    can freely move in the frame's x direction specified by the coordinate
    ``q``.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import (LinearPathway, Point,
    ...     ReferenceFrame)
    >>> from sympy.physics.vector import dynamicsymbols
    >>> N = ReferenceFrame('N')
    >>> q = dynamicsymbols('q')
    >>> force = Symbol('F')
    >>> pA, pB = Point('pA'), Point('pB')
    >>> pA.set_vel(N, 0)
    >>> pB.set_pos(pA, q * N.x)
    >>> pB.pos_from(pA)
    q(t)*N.x
    >>> linear_pathway = LinearPathway(pA, pB)
    >>> actuator = ForceActuator(force, linear_pathway)
    >>> actuator
    ForceActuator(F, LinearPathway(pA, pB))

    Parameters
    ==========

    force : Expr
        The scalar expression defining the (expansile) force that the actuator
        produces.
    pathway : PathwayBase
        The pathway that the actuator follows. This must be an instance of a
        concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.

    """

    def __init__(
        self,
        force: ExprType,
        pathway: PathwayBase,
    ) -> None:
        """Initializer for ``ForceActuator``.

        Parameters
        ==========

        force : Expr
            The scalar expression defining the (expansile) force that the
            actuator produces.
        pathway : PathwayBase
            The pathway that the actuator follows. This must be an instance of
            a concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.

        """
        self.force = force
        self.pathway = pathway

    @property
    def force(self) -> ExprType:
        """The magnitude of the force produced by the actuator."""
        return self._force

    @force.setter
    def force(self, force: ExprType) -> None:
        if hasattr(self, '_force'):
            msg = (
                f'Can\'t set attribute `force` to {repr(force)} as it is '
                f'immutable.'
            )
            raise AttributeError(msg)
        self._force = sympify(force, strict=True)

    @property
    def pathway(self) -> PathwayBase:
        """The ``Pathway`` defining the actuator's line of action."""
        return self._pathway

    @pathway.setter
    def pathway(self, pathway: PathwayBase) -> None:
        if hasattr(self, '_pathway'):
            msg = (
                f'Can\'t set attribute `pathway` to {repr(pathway)} as it is '
                f'immutable.'
            )
            raise AttributeError(msg)
        if not isinstance(pathway, PathwayBase):
            msg = (
                f'Value {repr(pathway)} passed to `pathway` was of type '
                f'{type(pathway)}, must be {PathwayBase}.'
            )
            raise TypeError(msg)
        self._pathway = pathway

    def to_loads(self) -> list[LoadBase]:
        """Loads required by the equations of motion method classes.

        Explanation
        ===========

        ``KanesMethod`` requires a list of ``Point``-``Vector`` tuples to be
        passed to the ``loads`` parameters of its ``kanes_equations`` method
        when constructing the equations of motion. This method acts as a
        utility to produce the correctly-structred pairs of points and vectors
        required so that these can be easily concatenated with other items in
        the list of loads and passed to ``KanesMethod.kanes_equations``. These
        loads are also in the correct form to also be passed to the other
        equations of motion method classes, e.g. ``LagrangesMethod``.

        Examples
        ========

        The below example shows how to generate the loads produced by a force
        actuator that follows a linear pathway. In this example we'll assume
        that the force actuator is being used to model a simple linear spring.
        First, create a linear pathway between two points separated by the
        coordinate ``q`` in the ``x`` direction of the global frame ``N``.

        >>> from sympy.physics.mechanics import (LinearPathway, Point,
        ...     ReferenceFrame)
        >>> from sympy.physics.vector import dynamicsymbols
        >>> q = dynamicsymbols('q')
        >>> N = ReferenceFrame('N')
        >>> pA, pB = Point('pA'), Point('pB')
        >>> pB.set_pos(pA, q * N.x)
        >>> pathway = LinearPathway(pA, pB)

        Now create a symbol ``k`` to describe the spring's stiffness and
        instantiate a force actuator that produces a (contractile) force
        proportional to both the spring's stiffness and the pathway's length.
        Note that actuator classes use the sign convention that expansile
        forces are positive, so for a spring to produce a contractile force the
        spring force needs to be calculated as the negative for the stiffness
        multiplied by the length.

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics._actuator import ForceActuator
        >>> stiffness = Symbol('k')
        >>> spring_force = -stiffness * pathway.length
        >>> spring = ForceActuator(spring_force, pathway)

        The forces produced by the spring can be generated in the list of loads
        form that ``KanesMethod`` (and other equations of motion methods)
        requires by calling the ``to_loads`` method.

        >>> spring.to_loads()
        [(pA, k*q(t)*N.x), (pB, - k*q(t)*N.x)]

        A simple linear damper can be modeled in a similar way. Create another
        symbol ``c`` to describe the dampers damping coefficient. This time
        instantiate a force actuator that produces a force proportional to both
        the damper's damping coefficient and the pathway's extension velocity.
        Note that the damping force is negative as it acts in the opposite
        direction to which the damper is changing in length.

        >>> damping_coefficient = Symbol('c')
        >>> damping_force = -damping_coefficient * pathway.extension_velocity
        >>> damper = ForceActuator(damping_force, pathway)

        Again, the forces produces by the damper can be generated by calling
        the ``to_loads`` method.

        >>> damper.to_loads()
        [(pA, c*Derivative(q(t), t)*N.x), (pB, - c*Derivative(q(t), t)*N.x)]

        """
        return self.pathway.compute_loads(self.force)

    def __repr__(self) -> str:
        """Representation of a ``ForceActuator``."""
        return f'{self.__class__.__name__}({self.force}, {self.pathway})'


class LinearSpring(ForceActuator):
    """A spring with its spring force as a linear function of its length.

    Explanation
    ===========

    Note that the "linear" in the name ``LinearSpring`` refers to the fact that
    the spring force is a linear function of the springs length. I.e. for a
    linear spring with stiffness ``k``, distance between its ends of ``x``, and
    an equilibrium length of ``0``, the spring force will be ``-k*x``, which is
    a linear function in ``x``. To create a spring that follows a linear, or
    straight, pathway between its two ends, a ``LinearPathway`` instance needs
    to be passed to the ``pathway`` parameter.

    Examples
    ========

    As the ``_actuator.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``LinearSpring`` must therefore be
    imported directly from the ``sympy.physics.mechanics._actuator`` module.

    >>> from sympy.physics.mechanics._actuator import LinearSpring

    To construct a linear spring, an expression (or symbol) must be supplied to
    represent the stiffness (spring constant) of the spring, alongside a
    pathway specifying its line of action. Let's also create a global reference
    frame and spatially fix one of the points in it while setting the other to
    be positioned such that it can freely move in the frame's x direction
    specified by the coordinate ``q``.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import (LinearPathway, Point,
    ...     ReferenceFrame)
    >>> from sympy.physics.vector import dynamicsymbols
    >>> N = ReferenceFrame('N')
    >>> q = dynamicsymbols('q')
    >>> stiffness = Symbol('k')
    >>> pA, pB = Point('pA'), Point('pB')
    >>> pA.set_vel(N, 0)
    >>> pB.set_pos(pA, q * N.x)
    >>> pB.pos_from(pA)
    q(t)*N.x
    >>> linear_pathway = LinearPathway(pA, pB)
    >>> spring = LinearSpring(stiffness, linear_pathway)
    >>> spring
    LinearSpring(k, LinearPathway(pA, pB))

    This spring will produce a force that is proportional to both its stiffness
    and the pathway's length. Note that this force is negative as SymPy's sign
    convention for actuators is that negative forces are contractile.

    >>> spring.force
    -k*sqrt(q(t)**2)

    To create a linear spring with a non-zero equilibrium length, an expression
    (or symbol) can be passed to the ``equilibrium_length`` parameter on
    construction on a ``LinearSpring`` instance. Let's create a symbol ``l``
    to denote a non-zero equilibrium length and create another linear spring.

    >>> l = Symbol('l')
    >>> spring = LinearSpring(stiffness, linear_pathway, equilibrium_length=l)
    >>> spring
    LinearSpring(k, LinearPathway(pA, pB), equilibrium_length=l)

    The spring force of this new spring is again proportional to both its
    stiffness and the pathway's length. However, the spring will not produce
    any force when ``q(t)`` equals ``l``. Note that the force will become
    expansile when ``q(t)`` is less than ``l``, as expected.

    >>> spring.force
    -k*(-l + sqrt(q(t)**2))

    Parameters
    ==========

    stiffness : Expr
        The spring constant.
    pathway : PathwayBase
        The pathway that the actuator follows. This must be an instance of a
        concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.
    equilibrium_length : Expr, optional
        The length at which the spring is in equilibrium, i.e. it produces no
        force. The default value is 0, i.e. the spring force is a linear
        function of the pathway's length with no constant offset.

    See Also
    ========

    ForceActuator: force-producing actuator (superclass of ``LinearSpring``).
    LinearPathway: straight-line pathway between a pair of points.

    """

    def __init__(
        self,
        stiffness: ExprType,
        pathway: PathwayBase,
        equilibrium_length: ExprType = S.Zero,
    ) -> None:
        """Initializer for ``LinearSpring``.

        Parameters
        ==========

        stiffness : Expr
            The spring constant.
        pathway : PathwayBase
            The pathway that the actuator follows. This must be an instance of
            a concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.
        equilibrium_length : Expr, optional
            The length at which the spring is in equilibrium, i.e. it produces
            no force. The default value is 0, i.e. the spring force is a linear
            function of the pathway's length with no constant offset.

        """
        self.stiffness = stiffness
        self.pathway = pathway
        self.equilibrium_length = equilibrium_length

    @property
    def force(self) -> ExprType:
        """The spring force produced by the linear spring."""
        return -self.stiffness * (self.pathway.length - self.equilibrium_length)

    @force.setter
    def force(self, force: Any) -> None:
        raise AttributeError('Can\'t set computed attribute `force`.')

    @property
    def stiffness(self) -> ExprType:
        """The spring constant for the linear spring."""
        return self._stiffness

    @stiffness.setter
    def stiffness(self, stiffness: ExprType):
        if hasattr(self, '_stiffness'):
            msg = (
                f'Can\'t set attribute `stiffness` to {repr(stiffness)} as it '
                f'is immutable.'
            )
            raise AttributeError(msg)
        self._stiffness = sympify(stiffness, strict=True)

    @property
    def equilibrium_length(self) -> ExprType:
        """The length of the spring at which it produces no force."""
        return self._equilibrium_length

    @equilibrium_length.setter
    def equilibrium_length(self, equilibrium_length: ExprType) -> None:
        if hasattr(self, '_equilibrium_length'):
            msg = (
                f'Can\'t set attribute `equilibrium_length` to '
                f'{repr(equilibrium_length)} as it is immutable.'
            )
            raise AttributeError(msg)
        self._equilibrium_length = sympify(equilibrium_length, strict=True)

    def __repr__(self) -> str:
        """Representation of a ``LinearSpring``."""
        string = f'{self.__class__.__name__}({self.stiffness}, {self.pathway}'
        if self.equilibrium_length == S.Zero:
            string += ')'
        else:
            string += f', equilibrium_length={self.equilibrium_length})'
        return string


class LinearDamper(ForceActuator):
    """A damper whose force is a linear function of its extension velocity.

    Explanation
    ===========

    Note that the "linear" in the name ``LinearDamper`` refers to the fact that
    the damping force is a linear function of the damper's rate of change in
    its length. I.e. for a linear damper with damping ``c`` and extension
    velocity ``v``, the damping force will be ``-c*v``, which is a linear
    function in ``v``. To create a damper that follows a linear, or straight,
    pathway between its two ends, a ``LinearPathway`` instance needs to be
    passed to the ``pathway`` parameter.

    Examples
    ========

    As the ``_actuator.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``LinearDamper`` must therefore be
    imported directly from the ``sympy.physics.mechanics._actuator`` module.

    >>> from sympy.physics.mechanics._actuator import LinearDamper

    To construct a linear damper, an expression (or symbol) must be supplied to
    represent the damping coefficient of the damper (we'll use the symbol
    ``c``), alongside a pathway specifying its line of action. Let's also
    create a global reference frame and spatially fix one of the points in it
    while setting the other to be positioned such that it can freely move in
    the frame's x direction specified by the coordinate ``q``. The velocity
    that the two points move away from one another can be specified by the
    coordinate ``u`` where ``u`` is the first time derivative of ``q``
    (i.e., ``u = Derivative(q(t), t)``).

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import (LinearPathway, Point,
    ...     ReferenceFrame)
    >>> from sympy.physics.vector import dynamicsymbols
    >>> N = ReferenceFrame('N')
    >>> q = dynamicsymbols('q')
    >>> damping = Symbol('c')
    >>> pA, pB = Point('pA'), Point('pB')
    >>> pA.set_vel(N, 0)
    >>> pB.set_pos(pA, q * N.x)
    >>> pB.pos_from(pA)
    q(t)*N.x
    >>> pB.vel(N)
    Derivative(q(t), t)*N.x
    >>> linear_pathway = LinearPathway(pA, pB)
    >>> damper = LinearDamper(damping, linear_pathway)
    >>> damper
    LinearDamper(c, LinearPathway(pA, pB))

    This damper will produce a force that is proportional to both its damping
    coefficient and the pathway's extension length. Note that this force is
    negative as SymPy's sign convention for actuators is that negative forces
    are contractile and the damping force of the damper will oppose the
    direction of length change.

    >>> damper.force
    -c*sqrt(q(t)**2)*Derivative(q(t), t)/q(t)

    Parameters
    ==========

    damping : Expr
        The damping constant.
    pathway : PathwayBase
        The pathway that the actuator follows. This must be an instance of a
        concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.

    See Also
    ========

    ForceActuator: force-producing actuator (superclass of ``LinearDamper``).
    LinearPathway: straight-line pathway between a pair of points.

    """

    def __init__(self, damping: ExprType, pathway: PathwayBase) -> None:
        """Initializer for ``LinearDamper``.

        Parameters
        ==========

        damping : Expr
            The damping constant.
        pathway : PathwayBase
            The pathway that the actuator follows. This must be an instance of
            a concrete subclass of ``PathwayBase``, e.g. ``LinearPathway``.

        """
        self.damping = damping
        self.pathway = pathway

    @property
    def force(self) -> ExprType:
        """The damping force produced by the linear damper."""
        return -self.damping * self.pathway.extension_velocity

    @force.setter
    def force(self, force: Any) -> None:
        raise AttributeError('Can\'t set computed attribute `force`.')

    @property
    def damping(self) -> ExprType:
        """The damping constant for the linear damper."""
        return self._damping

    @damping.setter
    def damping(self, damping: ExprType) -> None:
        if hasattr(self, '_damping'):
            msg = (
                f'Can\'t set attribute `damping` to {repr(damping)} as it is '
                f'immutable.'
            )
            raise AttributeError(msg)
        self._damping = sympify(damping, strict=True)

    def __repr__(self) -> str:
        """Representation of a ``LinearDamper``."""
        return f'{self.__class__.__name__}({self.damping}, {self.pathway})'


class TorqueActuator(ActuatorBase):
    """Torque-producing actuator.

    Explanation
    ===========

    A ``TorqueActuator`` is an actuator that produces a pair of equal and
    opposite torques on a pair of bodies.

    Examples
    ========

    As the ``_actuator.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``TorqueActuator`` must therefore be
    imported directly from the ``sympy.physics.mechanics._actuator`` module.

    >>> from sympy.physics.mechanics._actuator import TorqueActuator

    To construct a torque actuator, an expression (or symbol) must be supplied
    to represent the torque it can produce, alongside a vector specifying the
    axis about which the torque will act, and a pair of frames on which the
    torque will act.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, RigidBody
    >>> N = ReferenceFrame('N')
    >>> A = ReferenceFrame('A')
    >>> torque = Symbol('T')
    >>> axis = N.z
    >>> parent = RigidBody('parent', frame=N)
    >>> child = RigidBody('child', frame=A)
    >>> bodies = (child, parent)
    >>> actuator = TorqueActuator(torque, axis, *bodies)
    >>> actuator
    TorqueActuator(T, axis=N.z, target_frame=A, reaction_frame=N)

    Note that because torques actually act on frames, not bodies,
    ``TorqueActuator`` will extract the frame associated with a ``RigidBody``
    when one is passed instead of a ``ReferenceFrame``.

    Parameters
    ==========

    torque : Expr
        The scalar expression defining the torque that the actuator produces.
    axis : Vector
        The axis about which the actuator applies torques.
    target_frame : ReferenceFrame | RigidBody
        The primary frame on which the actuator will apply the torque.
    reaction_frame : ReferenceFrame | RigidBody | None
        The secondary frame on which the actuator will apply the torque. Note
        that the (equal and opposite) reaction torque is applied to this frame.

    """

    def __init__(
        self,
        torque: ExprType,
        axis: Vector,
        target_frame: ReferenceFrame | RigidBody,
        reaction_frame: ReferenceFrame | RigidBody | None = None,
    ) -> None:
        """Initializer for ``TorqueActuator``.

        Parameters
        ==========

        torque : Expr
            The scalar expression defining the torque that the actuator
            produces.
        axis : Vector
            The axis about which the actuator applies torques.
        target_frame : ReferenceFrame | RigidBody
            The primary frame on which the actuator will apply the torque.
        reaction_frame : ReferenceFrame | RigidBody | None
           The secondary frame on which the actuator will apply the torque.
           Note that the (equal and opposite) reaction torque is applied to
           this frame.

        """
        self.torque = torque
        self.axis = axis
        self.target_frame = target_frame  # type: ignore
        self.reaction_frame = reaction_frame  # type: ignore

    @classmethod
    def at_pin_joint(
        cls,
        torque: ExprType,
        pin_joint: PinJoint,
    ) -> TorqueActuator:
        """Alternate construtor to instantiate from a ``PinJoint`` instance.

        Examples
        ========

        To create a pin joint the ``PinJoint`` class requires a name, parent
        body, and child body to be passed to its constructor. It is also
        possible to control the joint axis using the ``joint_axis`` keyword
        argument. In this example let's use the parent body's reference frame's
        z-axis as the joint axis.

        >>> from sympy.physics.mechanics import (PinJoint, ReferenceFrame,
        ... RigidBody)
        >>> from sympy.physics.mechanics._actuator import TorqueActuator
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> parent = RigidBody('parent', frame=N)
        >>> child = RigidBody('child', frame=A)
        >>> pin_joint = PinJoint(
        ...     'pin',
        ...     parent,
        ...     child,
        ...     joint_axis=N.z,
        ... )

        Let's also create a symbol ``T`` that will represent the torque applied
        by the torque actuator.

        >>> from sympy import Symbol
        >>> torque = Symbol('T')

        To create the torque actuator from the ``torque`` and ``pin_joint``
        variables previously instantiated, these can be passed to the alternate
        constructor class method ``at_pin_joint`` of the ``TorqueActuator``
        class. It should be noted that a positive torque will cause a positive
        displacement of the joint coordinate or that the torque is applied on the
        child body with a reaction torque on the parent.

        >>> actuator = TorqueActuator.at_pin_joint(torque, pin_joint)
        >>> actuator
        TorqueActuator(T, axis=N.z, target_frame=A, reaction_frame=N)

        Parameters
        ==========

        torque : Expr
            The scalar expression defining the torque that the actuator
            produces.
        pin_joint : PinJoint
            The pin joint, and by association the parent and child bodies, on
            which the torque actuator will act. The pair of bodies acted upon
            by the torque actuator are the parent and child bodies of the pin
            joint, with the child acting as the reaction body. The pin joint's
            axis is used as the axis about which the torque actuator will apply
            its torque.

        """
        if not isinstance(pin_joint, PinJoint):
            msg = (
                f'Value {repr(pin_joint)} passed to `pin_joint` was of type '
                f'{type(pin_joint)}, must be {PinJoint}.'
            )
            raise TypeError(msg)
        return cls(
            torque,
            pin_joint.joint_axis,
            pin_joint.child_interframe,
            pin_joint.parent_interframe,
        )

    @property
    def torque(self) -> ExprType:
        """The magnitude of the torque produced by the actuator."""
        return self._torque

    @torque.setter
    def torque(self, torque: ExprType) -> None:
        if hasattr(self, '_torque'):
            msg = (
                f'Can\'t set attribute `torque` to {repr(torque)} as it is '
                f'immutable.'
            )
            raise AttributeError(msg)
        self._torque = sympify(torque, strict=True)

    @property
    def axis(self) -> Vector:
        """The axis about which the torque acts."""
        return self._axis

    @axis.setter
    def axis(self, axis: Vector) -> None:
        if hasattr(self, '_axis'):
            msg = (
                f'Can\'t set attribute `axis` to {repr(axis)} as it is '
                f'immutable.'
            )
            raise AttributeError(msg)
        if not isinstance(axis, Vector):
            msg = (
                f'Value {repr(axis)} passed to `axis` was of type '
                f'{type(axis)}, must be {Vector}.'
            )
            raise TypeError(msg)
        self._axis = axis

    @property
    def target_frame(self) -> ReferenceFrame:
        """The primary reference frames on which the torque will act."""
        return self._target_frame

    @target_frame.setter
    def target_frame(self, target_frame: ReferenceFrame) -> None:
        if hasattr(self, '_target_frame'):
            msg = (
                f'Can\'t set attribute `target_frame` to {repr(target_frame)} '
                f'as it is immutable.'
            )
            raise AttributeError(msg)
        if isinstance(target_frame, RigidBody):
            target_frame = target_frame.frame
        elif not isinstance(target_frame, ReferenceFrame):
            msg = (
                f'Value {repr(target_frame)} passed to `target_frame` was of '
                f'type {type(target_frame)}, must be {ReferenceFrame}.'
            )
            raise TypeError(msg)
        self._target_frame = target_frame

    @property
    def reaction_frame(self) -> ReferenceFrame | None:
        """The primary reference frames on which the torque will act."""
        return self._reaction_frame

    @reaction_frame.setter
    def reaction_frame(self, reaction_frame: ReferenceFrame | None) -> None:
        if hasattr(self, '_reaction_frame'):
            msg = (
                f'Can\'t set attribute `reaction_frame` to '
                f'{repr(reaction_frame)} as it is immutable.'
            )
            raise AttributeError(msg)
        if isinstance(reaction_frame, RigidBody):
            reaction_frame = reaction_frame.frame
        elif (
            not isinstance(reaction_frame, ReferenceFrame)
            and reaction_frame is not None
        ):
            msg = (
                f'Value {repr(reaction_frame)} passed to `reaction_frame` was '
                f'of type {type(reaction_frame)}, must be {ReferenceFrame}.'
            )
            raise TypeError(msg)
        self._reaction_frame = reaction_frame

    def to_loads(self) -> list[LoadBase]:
        """Loads required by the equations of motion method classes.

        Explanation
        ===========

        ``KanesMethod`` requires a list of ``Point``-``Vector`` tuples to be
        passed to the ``loads`` parameters of its ``kanes_equations`` method
        when constructing the equations of motion. This method acts as a
        utility to produce the correctly-structred pairs of points and vectors
        required so that these can be easily concatenated with other items in
        the list of loads and passed to ``KanesMethod.kanes_equations``. These
        loads are also in the correct form to also be passed to the other
        equations of motion method classes, e.g. ``LagrangesMethod``.

        Examples
        ========

        The below example shows how to generate the loads produced by a torque
        actuator that acts on a pair of bodies attached by a pin joint.

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import (PinJoint, ReferenceFrame,
        ... RigidBody)
        >>> from sympy.physics.mechanics._actuator import TorqueActuator
        >>> torque = Symbol('T')
        >>> N = ReferenceFrame('N')
        >>> A = ReferenceFrame('A')
        >>> parent = RigidBody('parent', frame=N)
        >>> child = RigidBody('child', frame=A)
        >>> pin_joint = PinJoint(
        ...     'pin',
        ...     parent,
        ...     child,
        ...     joint_axis=N.z,
        ... )
        >>> actuator = TorqueActuator.at_pin_joint(torque, pin_joint)

        The forces produces by the damper can be generated by calling the
        ``to_loads`` method.

        >>> actuator.to_loads()
        [(A, T*N.z), (N, - T*N.z)]

        Alternatively, if a torque actuator is created without a reaction frame
        then the loads returned by the ``to_loads`` method will contain just
        the single load acting on the target frame.

        >>> actuator = TorqueActuator(torque, N.z, N)
        >>> actuator.to_loads()
        [(N, T*N.z)]

        """
        loads: list[LoadBase] = [
            Torque(self.target_frame, self.torque * self.axis),
        ]
        if self.reaction_frame is not None:
            loads.append(Torque(self.reaction_frame, -self.torque * self.axis))
        return loads

    def __repr__(self) -> str:
        """Representation of a ``TorqueActuator``."""
        string = (
            f'{self.__class__.__name__}({self.torque}, axis={self.axis}, '
            f'target_frame={self.target_frame}'
        )
        if self.reaction_frame is not None:
            string += f', reaction_frame={self.reaction_frame})'
        else:
            string += ')'
        return string
