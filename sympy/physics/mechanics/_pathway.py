"""Implementations of pathways for use by actuators.

Notes
=====

This module is experimental and so is named with a leading underscore to
indicate that the API is not yet stabilized and could be subject to breaking
changes.

"""

from __future__ import annotations

from abc import ABC, abstractmethod

from sympy.core.backend import S, Symbol
from sympy.core.expr import Expr
from sympy.physics.mechanics import Point
from sympy.physics.vector import Vector, dynamicsymbols


__all__ = ['LinearPathway']


class PathwayBase(ABC):
    """Abstract base class for all pathway classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom pathway types through subclassing.

    """

    def __init__(
        self,
        origin: Point,
        insertion: Point,
    ) -> None:
        """Initializer for ``PathwayBase``."""
        self.origin = origin
        self.insertion = insertion

    @property
    def origin(self) -> Point:
        """The ``Point`` from which the pathway originates."""
        return self._origin

    @origin.setter
    def origin(self, origin: Point) -> None:
        if not isinstance(origin, Point):
            msg = (
                f'Value {repr(origin)} passed to `origin` was of type '
                f'{type(origin)}, must be {Point}.'
            )
            raise TypeError(msg)
        self._origin = origin

    @property
    def insertion(self) -> Point:
        """The ``Point`` into which the pathway inserts."""
        return self._insertion

    @insertion.setter
    def insertion(self, insertion: Point) -> None:
        if not isinstance(insertion, Point):
            msg = (
                f'Value {repr(insertion)} passed to `insertion` was of type '
                f'{type(insertion)}, must be {Point}.'
            )
            raise TypeError(msg)
        self._insertion = insertion

    @abstractmethod
    def _true_length(self) -> Expr:
        """Exact analytical expression for the pathway's length."""
        pass

    @abstractmethod
    def _true_shortening_velocity(self) -> Expr:
        """Exact analytical expression for the pathway's shortening velocity."""
        pass

    @property
    def length(self) -> Expr:
        """An expression representing the pathway's length."""
        return self._true_length()

    @property
    def shortening_velocity(self) -> Expr:
        """An expression representing the pathway's shortening velocity."""
        return self._true_shortening_velocity()

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}(origin={self.origin}, '
            f'insertion={self.insertion})'
        )


class LinearPathway(PathwayBase):
    """Linear pathway between origin and insertion.

    Explanation
    ===========

    A linear pathway forms a straight-line segment between two points and is
    the simplest pathway that can be formed. It will not interact with any
    other objects in the system, i.e. a ``LinearPathway`` will interest other
    objects to ensure that the path between its two ends (its origin and
    insertion) is the shortest possible.

    Examples
    ========

    As the ``_pathway.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace and must be imported directly from
    the ``sympy.physics.mechanics._pathway`` module.

    >>> from sympy.physics.mechanics._pathway import LinearPathway

    To construct a pathway, two points are required to be passed to the
    ``origin`` and ``insertion`` parameters.

    >>> from sympy.physics.mechanics import Point
    >>> origin = Point('pO')
    >>> insertion = Point('pI')
    >>> linear_pathway = LinearPathway(origin, insertion)
    >>> linear_pathway
    LinearPathway(origin=pO, insertion=pI)

    The pathway create above isn't very interesting without the positions and
    velocities of its origin and insertion points being described. Without this
    its not possible to describe how the pathway moves, i.e. its length or its
    shortening velocity.

    >>> from sympy.physics.mechanics import ReferenceFrame
    >>> from sympy.physics.vector import dynamicsymbols
    >>> N = ReferenceFrame('N')
    >>> q = dynamicsymbols('q')
    >>> insertion.set_pos(origin, q * N.x)
    >>> insertion.pos_from(origin)
    q(t)*N.x

    A pathway's length can be accessed via its ``length`` attribute.

    >>> linear_pathway.length
    sqrt(q(t)**2)

    Note how what appears to be an overly-complex expression is returned. This
    is actually required as it ensures that a pathway's length is always
    positive.

    A pathway's shortening velocity can be accessed similarly via its
    ``shortening_velocity`` attribute.

    >>> linear_pathway.shortening_velocity
    -q(t)*Derivative(q(t), t)/sqrt(q(t)**2)

    Parameters
    ==========

    origin : Point
        The ``Point`` from which the pathway originates.
    insertion : Point
        The ``Point`` into which the pathway inserts.

    """

    def __init__(
        self,
        origin: Point,
        insertion: Point,
    ) -> None:
        """Initializer for ``LinearPathway``.

        Parameters
        ==========

        origin : Point
            The ``Point`` from which the pathway originates.
        insertion : Point
            The ``Point`` into which the pathway inserts.

        """
        super().__init__(
            origin=origin,
            insertion=insertion,
        )

    def _true_length(self) -> Expr:
        """Exact analytical expression for the pathway's length."""
        length = self.insertion.pos_from(self.origin).magnitude()
        return length

    def _true_shortening_velocity(self) -> Expr:
        """Exact analytical expression for the pathway's shortening velocity."""
        relative_position = self.insertion.pos_from(self.origin)
        if not relative_position:
            return S.Zero
        t = dynamicsymbols._t  # type: ignore
        # A reference frame is needed to differentiate ``relative_position`` to
        # ``relative_velocity`` so choose the first ``ReferenceFrame`` that
        # ``relative_position`` is defined using.
        frame = relative_position.args[0][1]
        relative_velocity = relative_position.diff(t, frame)
        shortening_velocity = -relative_velocity.dot(relative_position.normalize())
        return shortening_velocity

    def forces(self, force: Symbol) -> list[tuple[Point, Vector]]:
        """Forces list required by ``KanesMethod``.

        Explanation
        ===========

        ``KanesMethod`` requires a ``list[tuple[Point, Vector]]`` to be passed
        to its ``forcelist`` parameter on creation of an instance. This method
        acts as a utility to produce the correctly-structred pairs of points
        and vectors required so that these can be easily concatenated with
        other items in the force list and passed to ``KanesMethod``.

        Examples
        ========

        The below example shows how to generate the forces produced in a linear
        actuator that produces a contractile force ``F``. First, create a
        linear actuator between two points separated by the coordinate ``q``
        in the ``x`` direction of the global frame ``N``.

        >>> from sympy.physics.mechanics import Point, ReferenceFrame
        >>> from sympy.physics.mechanics._pathway import LinearPathway
        >>> from sympy.physics.vector import dynamicsymbols
        >>> q = dynamicsymbols('q')
        >>> N = ReferenceFrame('N')
        >>> origin = Point('pO')
        >>> insertion = Point('pI')
        >>> insertion.set_pos(origin, q * N.x)
        >>> linear_pathway = LinearPathway(origin, insertion)

        Now create a symbol ``F`` to describe the magnitude of the
        (contractile) for that will be produced along the pathway. The list of
        forces that ``KanesMethod`` requires can be produced by calling the
        pathway's ``forces`` method with ``F`` passed as the only argument.

        >>> from sympy import Symbol
        >>> F = Symbol('F')
        >>> linear_pathway.forces(F)
        [(pO, F*q(t)/sqrt(q(t)**2)*N.x), (pI, - F*q(t)/sqrt(q(t)**2)*N.x)]

        Parameters
        ==========

        force : Symbol
            The force produced by the actuator. It is assumed that this
            ``Symbol`` represents a contractile force.

        """
        forces = [
            (
                self.origin,
                force * self.insertion.pos_from(self.origin) / self.length,
            ),
            (   self.insertion,
                force * self.origin.pos_from(self.insertion) / self.length,
            ),
        ]
        return forces
