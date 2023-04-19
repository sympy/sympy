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

from sympy.core.backend import Integer, Mul, acos, pi, sqrt, tan
from sympy.functions.elementary.trigonometric import atan2
from sympy.polys.polytools import cancel
from sympy.simplify.simplify import posify, trigsimp

if TYPE_CHECKING:
    from sympy.core.backend import Symbol
    from sympy.core.expr import Expr
    from sympy.physics.mechanics import Point, Vector


__all__ = ['Sphere']


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

    @abstractmethod
    def geodesic_length(self, point_1: Point, point_2: Point) -> Expr:
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
        The radius of the sphere.
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
    def radius(self, radius) -> None:
        self._radius = radius

    @property
    def point(self) -> Point:
        """A point through which the cylinder's axis passes."""
        return self._point

    @point.setter
    def point(self, point) -> None:
        self._point = point

    def _point_is_on_surface(self, point: Point) -> bool:
        """Determine if a point is on the sphere's surface.

        Parameters
        ==========

        point : Point
            The point for which it's to be ascertained if it's on the sphere's
            surface or not.

        """
        point_radius = _cancel_sqrt_of_squared_positive(
            point.pos_from(self.point).magnitude()
        )
        return point_radius == self.radius

    def geodesic_length(self, point_1: Point, point_2: Point) -> Expr:
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

        The geodesic length, which is in this case a quarter of the sphere's
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
        central_angle = acos(cancel(point_2_vector.dot(point_1_vector)))
        geodesic_length = self.radius * central_angle
        return geodesic_length

    def __repr__(self) -> str:
        """Representation of a ``Sphere``."""
        return (
            f'{self.__class__.__name__}(radius={self.radius}, '
            f'point={self.point})'
        )


def _cancel_sqrt_of_squared_positive(expr) -> Expr:
    """Cancel ``sqrt(x**2)`` to ``x`` in expressions for positive ``x``."""
    expr, reps = posify(expr)
    return expr.subs(reps)
