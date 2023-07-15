"""Implementations of geometry objects for use by wrapping pathways.

Notes
=====

This module is experimental and so is named with a leading underscore to
indicate that the API is not yet stabilized and could be subject to breaking
changes.

"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from sympy.core.backend import Integer, acos, pi, sqrt, sympify, tan
from sympy.core.relational import Eq
from sympy.functions.elementary.trigonometric import atan2
from sympy.polys.polytools import cancel
from sympy.physics.vector import Vector, dot
from sympy.simplify.simplify import trigsimp

if TYPE_CHECKING:
    from sympy.core.backend import USE_SYMENGINE, Symbol
    from sympy.physics.mechanics import Point

    if USE_SYMENGINE:
        from sympy.core.backend import Basic as ExprType
    else:
        from sympy.core.expr import Expr as ExprType


__all__ = [
    'Cylinder',
    'Sphere',
]


class GeometryBase(ABC):
    """Abstract base class for all geometry classes to inherit from.

    Notes
    =====

    Instances of this class cannot be directly instantiated by users. However,
    it can be used to created custom geometry types through subclassing.

    """

    @abstractmethod
    def _point_is_on_surface(self, point: Point) -> bool:
        """Determine if a point is on the geometry's surface.

        Parameters
        ==========
        point : Point
            The point for which it's to be ascertained if it's on the
            geometry's surface or not.

        """
        pass

    @abstractmethod
    def geodesic_length(self, point_1: Point, point_2: Point) -> ExprType:
        """The shortest distance between two points on a geometry's surface.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic length should be calculated.
        point_2 : Point
            The point to which the geodesic length should be calculated.

        """
        pass

    def __repr__(self) -> str:
        """Default representation of a geometry object."""
        return f'{self.__class__.__name__}()'


class Sphere(GeometryBase):
    """A solid spherical object.

    Explanation
    ===========

    A wrapping geometry that allows for circular arcs to be defined between
    pairs of points. These paths are always geodetic (the shortest possible).

    Examples
    ========

    As the ``_geometry.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``Sphere`` must therefore be
    imported directly from the ``sympy.physics.mechanics._geometry`` module.

    >>> from sympy.physics.mechanics._geometry import Sphere

    To create a ``Sphere`` instance, a ``Symbol`` denoting its radius and
    ``Point`` at which its center will be located are needed:

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Point
    >>> r = Symbol('r')
    >>> pO = Point('pO')


    A sphere with radius ``r`` centered on ``pO`` can be instantiated with:

    >>> Sphere(r, pO)
    Sphere(radius=r, point=pO)

    Parameters
    ==========

    radius : Symbol
        The radius of the sphere. This symbol must represent a value that is
        positive and constant, i.e. it cannot be a dynamic symbol.
    point : Point
        A point at which the sphere is centered.

    See Also
    ========

    Cylinder: Cylindrical geometry where the wrapping direction can be defined.

    """

    def __init__(self, radius: Symbol, point: Point) -> None:
        """Initializer for ``Sphere``.

        Parameters
        ==========

        radius : Symbol
            The radius of the sphere.
        point : Point
            A point through which the sphere's axis passes.

        """
        self.radius = radius
        self.point = point

    @property
    def radius(self) -> Symbol:
        """The radius of the sphere."""
        return self._radius

    @radius.setter
    def radius(self, radius: Symbol) -> None:
        self._radius = radius

    @property
    def point(self) -> Point:
        """A point through which the cylinder's axis passes."""
        return self._point

    @point.setter
    def point(self, point: Point) -> None:
        self._point = point

    def _point_is_on_surface(self, point: Point) -> bool:
        """Determine if a point is on the sphere's surface.

        Parameters
        ==========

        point : Point
            The point for which it's to be ascertained if it's on the sphere's
            surface or not.

        """
        point_vector = point.pos_from(self.point)
        if isinstance(point_vector, Vector):
            point_radius = dot(point_vector, point_vector)
        else:
            point_radius = point_vector**2
        return Eq(point_radius, self.radius**2) == True

    def geodesic_length(self, point_1: Point, point_2: Point) -> ExprType:
        r"""The shortest distance between two points on a geometry's surface.

        Explanation
        ===========

        The geodesic length, i.e. the shortest arc along the surface of a
        sphere, connecting two points can be calculated using the formula:

        .. math::
            l = \acos{\mathbf{v}_1 \dot \mathbf{v}_2}

        where $mathbf{v}_1$ and $mathbf{v}_2$ are the unit vectors from the
        sphere's center to the first and second points on the sphere's surface
        respectively.

        Note that the actual path that the geodesic will take is undefined when
        the two points are directly opposite one another.

        Examples
        ========

        A geodesic length can only be calculated between two points on the
        sphere's surface. Firstly, a ``Sphere`` instance must be created along
        with two points that will lie on its surface:

        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Point, ReferenceFrame
        >>> from sympy.physics.mechanics._geometry import Sphere
        >>> N = ReferenceFrame('N')
        >>> r = Symbol('r')
        >>> pO = Point('pO')
        >>> pO.set_vel(N, 0)
        >>> sphere = Sphere(r, pO)
        >>> p1 = Point('p1')
        >>> p2 = Point('p2')

        Let's assume that ``p1`` lies at a distance of ``r`` in the ``N.x``
        direction from ``pO`` and that ``p2`` is located on the sphere's
        surface in the ``N.y + N.z`` direction from ``pO``. These positions can
        be set with:

        >>> p1.set_pos(pO, r * N.x)
        >>> p1.pos_from(pO)
        r*N.x
        >>> p2.set_pos(pO, r * (N.y + N.z).normalize())
        >>> p2.pos_from(pO)
        sqrt(2)*r/2*N.y + sqrt(2)*r/2*N.z

        The geodesic length, which is in this case is a quarter of the sphere's
        circumference, can be calculated using the ``geodesic_length`` method:

        >>> sphere.geodesic_length(p1, p2)
        pi*r/2

        If the ``geodesic_length`` method is passed an argument ``Point`` that
        doesn't lie on the sphere's surface then a ``ValueError`` is raised
        because it's not possible to calculate a value in this case.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic length should be calculated.
        point_2 : Point
            The point to which the geodesic length should be calculated.

        """
        for point in (point_1, point_2):
            if not self._point_is_on_surface(point):
                msg = (
                    f'Geodesic length cannot be calculated as point {point} '
                    f'with radius {point.pos_from(self.point).magnitude()} '
                    f'from the sphere\'s center {self.point} does not lie on '
                    f'the surface of {self} with radius {self.radius}.'
                )
                raise ValueError(msg)
        point_1_vector = point_1.pos_from(self.point).normalize()
        point_2_vector = point_2.pos_from(self.point).normalize()
        central_angle = acos(point_2_vector.dot(point_1_vector))
        geodesic_length = self.radius * central_angle
        return geodesic_length

    def __repr__(self) -> str:
        """Representation of a ``Sphere``."""
        return (
            f'{self.__class__.__name__}(radius={self.radius}, '
            f'point={self.point})'
        )


class Cylinder(GeometryBase):
    """A solid (infinite) cylindrical object.

    Explanation
    ===========

    A wrapping geometry that allows for circular arcs to be defined between
    pairs of points. These paths are always geodetic (the shortest possible) in
    the sense that they will be a straight line on the unwrapped cylinder's
    surface. However, it is also possible for a direction to be specified, i.e.
    paths can be influenced such that they either wrap along the shortest side
    or the longest side of the cylinder. To define these directions, rotations
    are in the positive direction following the right-hand rule.

    Examples
    ========

    As the ``_geometry.py`` module is experimental, it is not yet part of the
    ``sympy.physics.mechanics`` namespace. ``Cylinder`` must therefore be
    imported directly from the ``sympy.physics.mechanics._geometry`` module.

    >>> from sympy.physics.mechanics._geometry import Cylinder

    To create a ``Cylinder`` instance, a ``Symbol`` denoting its radius, a
    ``Vector`` defining its axis, and a ``Point`` through which its axis passes
    are needed:

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Point, ReferenceFrame
    >>> N = ReferenceFrame('N')
    >>> r = Symbol('r')
    >>> pO = Point('pO')
    >>> ax = N.x

    A cylinder with radius ``r``, and axis parallel to ``N.x`` passing through
    ``pO`` can be instantiated with:

    >>> Cylinder(r, pO, ax)
    Cylinder(radius=r, point=pO, axis=N.x)

    Parameters
    ==========

    radius : Symbol
        The radius of the cylinder.
    point : Point
        A point through which the cylinder's axis passes.
    axis : Vector
        The axis along which the cylinder is aligned.

    See Also
    ========

    Sphere: Spherical geometry where the wrapping direction is always geodetic.

    """

    def __init__(self, radius: Symbol, point: Point, axis: Vector) -> None:
        """Initializer for ``Cylinder``.

        Parameters
        ==========

        radius : Symbol
            The radius of the cylinder. This symbol must represent a value that
            is positive and constant, i.e. it cannot be a dynamic symbol.
        point : Point
            A point through which the cylinder's axis passes.
        axis : Vector
            The axis along which the cylinder is aligned.

        """
        self.radius = radius
        self.point = point
        self.axis = axis

    @property
    def radius(self) -> Symbol:
        """The radius of the cylinder."""
        return self._radius

    @radius.setter
    def radius(self, radius: Symbol) -> None:
        self._radius = radius

    @property
    def point(self) -> Point:
        """A point through which the cylinder's axis passes."""
        return self._point

    @point.setter
    def point(self, point: Point) -> None:
        self._point = point

    @property
    def axis(self) -> Vector:
        """The axis along which the cylinder is aligned."""
        return self._axis

    @axis.setter
    def axis(self, axis: Vector) -> None:
        self._axis = axis

    def _point_is_on_surface(self, point: Point) -> bool:
        """Determine if a point is on the cylinder's surface.

        Parameters
        ==========

        point : Point
            The point for which it's to be ascertained if it's on the
            cylinder's surface or not.

        """
        relative_position = point.pos_from(self.point)
        parallel = relative_position.dot(self.axis) * self.axis
        point_vector = relative_position - parallel
        if isinstance(point_vector, Vector):
            point_radius = dot(point_vector, point_vector)
        else:
            point_radius = point_vector**2
        return Eq(trigsimp(point_radius), self.radius**2) == True

    def geodesic_length(self, point_1: Point, point_2: Point) -> ExprType:
        r"""The shortest distance between two points on a geometry's surface.

        Explanation
        ===========

        The geodesic length, i.e. the shortest arc along the surface of a
        cylinder, connecting two points can be calculated using the formula:

        .. math::
            l = \acos{\mathbf{v}_1 \dot \mathbf{v}_2}

        where $mathbf{v}_1$ and $mathbf{v}_2$ are the unit vectors from the
        cylinder's center to the first and second points on the cylinder's
        surface respectively.

        Examples
        ========

        A geodesic length can only be calculated between two points on the
        cylinder's surface. Firstly, a ``Cylinder`` instance must be created
        along with two points that will lie on its surface:

        >>> from sympy import Symbol, cos, sin
        >>> from sympy.physics.mechanics import (Point, ReferenceFrame,
        ... dynamicsymbols)
        >>> from sympy.physics.mechanics._geometry import Cylinder
        >>> N = ReferenceFrame('N')
        >>> r = Symbol('r')
        >>> pO = Point('pO')
        >>> pO.set_vel(N, 0)
        >>> cylinder = Cylinder(r, pO, N.x)
        >>> p1 = Point('p1')
        >>> p2 = Point('p2')

        Let's assume that ``p1`` is located at ``N.x + r*N.y`` relative to
        ``pO`` and that ``p2`` is located at ``r*(cos(q)*N.y + sin(q)*N.z)``
        relative to ``pO``, where ``q(t)`` is a generalized coordinate
        specifying the angle rotated around the ``N.x`` axis according to the
        right-hand rule where ``N.y`` is zero. These positions can be set with:

        >>> q = dynamicsymbols('q')
        >>> p1.set_pos(pO, N.x + r*N.y)
        >>> p1.pos_from(pO)
        N.x + r*N.y
        >>> p2.set_pos(pO, r * (cos(q)*N.y + sin(q)*N.z).normalize())
        >>> p2.pos_from(pO).simplify()
        r*cos(q(t))*N.y + r*sin(q(t))*N.z

        The geodesic length, which is in this case a is the hypotenuse of a
        right triangle where the other two side lengths are ``1`` (parallel to
        the cylinder's axis) and ``r*q(t)`` (parallel to the cylinder's cross
        section), can be calculated using the ``geodesic_length`` method:

        >>> cylinder.geodesic_length(p1, p2).simplify()
        sqrt(r**2*q(t)**2 + 1)

        If the ``geodesic_length`` method is passed an argument ``Point`` that
        doesn't lie on the sphere's surface then a ``ValueError`` is raised
        because it's not possible to calculate a value in this case.

        Parameters
        ==========

        point_1 : Point
            The point from which the geodesic length should be calculated.
        point_2 : Point
            The point to which the geodesic length should be calculated.

        """
        for point in (point_1, point_2):
            if not self._point_is_on_surface(point):
                msg = (
                    f'Geodesic length cannot be calculated as point {point} '
                    f'with radius {point.pos_from(self.point).magnitude()} '
                    f'from the cylinder\'s center {self.point} does not lie on '
                    f'the surface of {self} with radius {self.radius} and axis '
                    f'{self.axis}.'
                )
                raise ValueError(msg)

        relative_position = point_2.pos_from(point_1)
        parallel_length = relative_position.dot(self.axis)

        point_1_relative_position = point_1.pos_from(self.point)
        point_1_perpendicular_vector = (
            point_1_relative_position
            - point_1_relative_position.dot(self.axis) * self.axis
        ).normalize()

        point_2_relative_position = point_2.pos_from(self.point)
        point_2_perpendicular_vector = (
            point_2_relative_position
            - point_2_relative_position.dot(self.axis) * self.axis
        ).normalize()

        central_angle = _directional_atan(
            cancel(point_1_perpendicular_vector
                .cross(point_2_perpendicular_vector)
                .dot(self.axis)),
            cancel(point_1_perpendicular_vector.dot(point_2_perpendicular_vector)),
        )

        planar_arc_length = self.radius * central_angle
        geodesic_length = sqrt(parallel_length**2 + planar_arc_length**2)
        return geodesic_length

    def __repr__(self) -> str:
        """Representation of a ``Cylinder``."""
        return (
            f'{self.__class__.__name__}(radius={self.radius}, '
            f'point={self.point}, axis={self.axis})'
        )


def _directional_atan(numerator: ExprType, denominator: ExprType) -> ExprType:
    """Compute atan in a directional sense as required for geodesics.

    Explanation
    ===========

    To be able to control the direction of the geodesic length along the
    surface of a cylinder a dedicated arctangent function is needed that
    properly handles the directionality of different case. This function
    ensures that the central angle is always positive but shifting the case
    where ``atan2`` would return a negative angle to be centered around
    ``2*pi``.

    Notes
    =====

    This function only handles very specific cases, i.e. the ones that are
    expected to be encountered when calculating symbolic geodesics on uniformly
    curved surfaces. As such, ``NotImplemented`` errors can be raised in many
    cases. This function is named with a leader underscore to indicate that it
    only aims to provide very specific functionality within the private scope
    of this module.

    """

    if numerator.is_number and denominator.is_number:
        angle = atan2(numerator, denominator)
        if angle < 0:
            angle += 2 * pi
    elif numerator.is_number:
        msg = (
            f'Cannot compute a directional atan when the numerator {numerator} '
            f'is numeric and the denominator {denominator} is symbolic.'
        )
        raise NotImplementedError(msg)
    elif denominator.is_number:
        msg = (
            f'Cannot compute a directional atan when the numerator {numerator} '
            f'is symbolic and the denominator {denominator} is numeric.'
        )
        raise NotImplementedError(msg)
    else:
        ratio = sympify(trigsimp(numerator / denominator))
        if isinstance(ratio, tan):
            angle = ratio.args[0]
        elif (
            ratio.is_Mul
            and ratio.args[0] == Integer(-1)
            and isinstance(ratio.args[1], tan)
        ):
            angle = 2 * pi - ratio.args[1].args[0]
        else:
            msg = f'Cannot compute a directional atan for the value {ratio}.'
            raise NotImplementedError(msg)

    return angle
